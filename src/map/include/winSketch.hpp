/**
 * @file    winSketch.hpp
 * @brief   routines to index the reference 
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef WIN_SKETCH_HPP 
#define WIN_SKETCH_HPP

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <cassert>
#include <zlib.h>  

//Own includes
#include "map/include/base_types.hpp"
#include "map/include/map_parameters.hpp"
#include "map/include/commonFunc.hpp"
#include "map/include/ThreadPool.hpp"

//External includes
#include "common/csv.h"
#include "common/kseq.h"
#include "common/murmur3.h"
#include "common/prettyprint.hpp"
#include "common/sparsehash/dense_hash_map"

KSEQ_INIT(gzFile, gzread)

namespace skch
{
  /**
   * @class     skch::Sketch
   * @brief     sketches and indexes the reference (subject sequence)
   * @details  
   *            1.  Minimizers are computed in streaming fashion
   *                Computing minimizers is using double ended queue which gives
   *                O(reference size) complexity
   *                Algorithm described here:
   *                https://people.cs.uct.ac.za/~ksmith/articles/sliding_window_minimum.html
   *
   *            2.  Index hashes into appropriate format to enable fast search at L1 mapping stage
   */
  class Sketch
    {
      //private members
    
      //algorithm parameters
      const skch::Parameters &param;

      //Ignore top % most frequent minimizers while lookups
      const float percentageThreshold = 0.001;

      //Minimizers that occur this or more times will be ignored (computed based on percentageThreshold)
      int freqThreshold = std::numeric_limits<int>::max();

      //Set of frequent seeds to be ignored
      std::unordered_set<hash_t> frequentSeeds;

      //Make the default constructor private, non-accessible
      Sketch();

      public:

      typedef std::vector< MashimizerInfo > MI_Type;
      using MIIter_t = MI_Type::const_iterator;

      //Keep sequence length, name that appear in the sequence (for printing the mappings later)
      std::vector< ContigInfo > metadata;

      /*
       * Keep the information of what sequences come from what file#
       * Example [a, b, c] implies 
       *  file 0 contains 0 .. a-1 sequences
       *  file 1 contains a .. b-1 
       *  file 2 contains b .. c-1
       */
      std::vector< int > sequencesByFileInfo;

      //Index for fast seed lookup (unordered_map)
      /*
       * [minimizer #1] -> [pos1, pos2, pos3 ...]
       * [minimizer #2] -> [pos1, pos2...]
       * ...
       */
      //using MI_Map_t = google::dense_hash_map< MinimizerMapKeyType, MinimizerMapValueType >;
      using MI_Map_t = std::unordered_map< MinimizerMapKeyType, MinimizerMapValueType >;
      MI_Map_t minimizerPosLookupIndex;

      /**
       * Keep list of minimizers, sequence# , their position within seq , here while parsing sequence 
       * Note : position is local within each contig
       * Hashes saved here are non-unique, ordered as they appear in the reference
       */
      MI_Type mashimizerIndex;

      private:


      //Frequency histogram of minimizers
      //[... ,x -> y, ...] implies y number of minimizers occur x times
      std::map<int, int> minimizerFreqHistogram;

      public:

      /**
       * @brief   constructor
       *          also builds, indexes the minimizer table
       */
      Sketch(const skch::Parameters &p) 
        :
          param(p) {
            if (p.indexFileName != "") {
              this->build(p.indexFileName);
            } else {
              this->build();
            }
            this->index();
            this->computeFreqHist();
            this->computeFreqSeedSet();
            if (p.outIndex != "") {
              this->write(p.outIndex);
            }
            // TODO should write these out and save instead of drop every time
            this->dropFreqSeedSet();
          }

      private:

      /**
       * @brief     build the sketch table
       * @details   compute and save minimizers from the reference sequence(s)
       *            assuming a fixed window size
       */
      void build(std::string indexFileName)
      {
        this->read(indexFileName);
        seqno_t seqCounter = 0;
        for(const auto &fileName : param.refSequences)
        {
          //Open the file using kseq
          gzFile fp = gzopen(fileName.c_str(), "r");
          kseq_t *seq = kseq_init(fp);


          //size of sequence
          offset_t len;

          while ((len = kseq_read(seq)) >= 0) 
          {
            //Save the sequence name
            metadata.push_back( ContigInfo{seq->name.s, (offset_t)seq->seq.l} );

            //Is the sequence too short?
            if(len < param.windowSize || len < param.kmerSize)
            {
#ifdef DEBUG
              std::cout << "WARNING, skch::Sketch::build, found an unusually short sequence relative to kmer and window size" << std::endl;
#endif
              seqCounter++;
              continue;  
            }
            else
            seqCounter++;
          }

          sequencesByFileInfo.push_back(seqCounter);

          kseq_destroy(seq);  
          gzclose(fp); //close the file handler 
        }

      }

      /**
       * @brief     build the sketch table
       * @details   compute and save minimizers from the reference sequence(s)
       *            assuming a fixed window size
       */
      void build()
      {
        //sequence counter while parsing file
        seqno_t seqCounter = 0;

        //Create the thread pool 
        ThreadPool<InputSeqContainer, MI_Type> threadPool( [this](InputSeqContainer* e) {return buildHelper(e);}, param.threads);

        for(const auto &fileName : param.refSequences)
        {

#ifdef DEBUG
        std::cout << "INFO, skch::Sketch::build, building minimizer index for " << fileName << std::endl;
#endif

          //Open the file using kseq
          gzFile fp = gzopen(fileName.c_str(), "r");
          kseq_t *seq = kseq_init(fp);


          //size of sequence
          offset_t len;

          while ((len = kseq_read(seq)) >= 0) 
          {
            //Save the sequence name
            metadata.push_back( ContigInfo{seq->name.s, (offset_t)seq->seq.l} );

            //Is the sequence too short?
            if(len < param.windowSize || len < param.kmerSize)
            {
#ifdef DEBUG
              std::cout << "WARNING, skch::Sketch::build, found an unusually short sequence relative to kmer and window size" << std::endl;
#endif
              seqCounter++;
              continue;  
            }
            else
            {
              threadPool.runWhenThreadAvailable(new InputSeqContainer(seq->seq.s, seq->name.s, len, seqCounter));

              //Collect output if available
              while ( threadPool.outputAvailable() )
                this->buildHandleThreadOutput(threadPool.popOutputWhenAvailable());
            }

            seqCounter++;
          }

          sequencesByFileInfo.push_back(seqCounter);

          kseq_destroy(seq);  
          gzclose(fp); //close the file handler 
        }

        //Collect remaining output objects
        while ( threadPool.running() )
          this->buildHandleThreadOutput(threadPool.popOutputWhenAvailable());

        std::cout << "INFO, skch::Sketch::build, minimizers picked from reference = " << mashimizerIndex.size() << std::endl;

      }

      /**
       * @brief               function to compute minimizers given input sequence object
       * @details             this function is run in parallel by multiple threads
       * @param[in]   input   input read details
       * @return              output object containing the mappings
       */
      MI_Type* buildHelper(InputSeqContainer *input)
      {
        MI_Type* thread_output = new MI_Type();

        //Compute minimizers in reference sequence
        skch::CommonFunc::addMashimizers(
                *thread_output, 
                &(input->seq[0u]), 
                input->len, 
                param.kmerSize, 
                param.segLength, 
                param.alphabetSize, 
                param.sketchSize,
                input->seqCounter);

        return thread_output;
      }

      /**
       * @brief                 routine to handle thread's local minimizer index
       * @param[in] output      thread local minimizer output
       */
      void buildHandleThreadOutput(MI_Type* output)
      {
        this->mashimizerIndex.insert(this->mashimizerIndex.end(), output->begin(), output->end());
        delete output;
      }

      /**
       * @brief   build the index for fast lookups using minimizer table
       */
      void index()
      {
        //Parse all the minimizers and push into the map
        //minimizerPosLookupIndex.set_empty_key(0);

        for(auto &e : mashimizerIndex)
        {
          minimizerPosLookupIndex[e.hash].push_back(e);
        }

        std::cout << "INFO, skch::Sketch::index, unique minimizers = " << minimizerPosLookupIndex.size() << std::endl;
      }

      void write(std::string indexFileName) {
        std::ofstream outStream;
        outStream.open(indexFileName);
        outStream << "seqId" << "\t" << "strand" << "\t" << "start" << "\t" << "end" << "\t" << "hash" << std::endl;
        for (auto& mi : mashimizerIndex) {
          outStream << mi.seqId << "\t" << std::to_string(mi.strand) << "\t" << mi.wpos << "\t" << mi.wpos_end << "\t" << mi.hash << std::endl;
        }
        outStream.close(); 
      }

      void read(std::string indexFileName) {
        io::CSVReader<5, io::trim_chars<' '>, io::no_quote_escape<'\t'>> inReader(indexFileName);
        inReader.read_header(io::ignore_missing_column, "seqId", "strand", "start", "end", "hash");
        hash_t hash;
        offset_t start, end;
        strand_t strand;
        seqno_t seqId;
        while (inReader.read_row(seqId, strand, start, end, hash))
        {
          mashimizerIndex.push_back(MashimizerInfo {hash, seqId, start, end, strand});
        }
      }

      /**
       * @brief   report the frequency histogram of minimizers using position lookup index
       *          and compute which high frequency minimizers to ignore
       */
      void computeFreqHist()
      {

        //1. Compute histogram

        for(auto &e : this->minimizerPosLookupIndex)
          this->minimizerFreqHistogram[e.second.size()] += 1;

        std::cout << "INFO, skch::Sketch::computeFreqHist, Frequency histogram of minimizers = " <<  *this->minimizerFreqHistogram.begin() <<  " ... " << *this->minimizerFreqHistogram.rbegin() << std::endl;

        //2. Compute frequency threshold to ignore most frequent minimizers

        int64_t totalUniqueMinimizers = this->minimizerPosLookupIndex.size();
        int64_t minimizerToIgnore = totalUniqueMinimizers * percentageThreshold / 100;

        int64_t sum = 0;

        //Iterate from highest frequent minimizers
        for(auto it = this->minimizerFreqHistogram.rbegin(); it != this->minimizerFreqHistogram.rend(); it++)
        {
          sum += it->second; //add frequency
          if(sum < minimizerToIgnore)
          {
            this->freqThreshold = it->first;
            //continue
          }
          else if(sum == minimizerToIgnore)
          {
            this->freqThreshold = it->first;
            break;
          }
          else
          {
            break;
          }
        }

        if(this->freqThreshold != std::numeric_limits<int>::max())
          std::cout << "INFO, skch::Sketch::computeFreqHist, With threshold " << this->percentageThreshold << "%, ignore minimizers occurring >= " << this->freqThreshold << " times during lookup." << std::endl;
        else
          std::cout << "INFO, skch::Sketch::computeFreqHist, With threshold " << this->percentageThreshold << "%, consider all minimizers during lookup." << std::endl;

      }

      /**
       * @brief     Construct a set of all hashes that appear above the threshold
       */
      void computeFreqSeedSet()
      {
        for(auto &e : this->minimizerPosLookupIndex) {
          if (e.second.size() >= this->freqThreshold) {
            this->frequentSeeds.insert(e.first);
          }
        }
      }

      /**
       * @brief     Remove frequent seeds from set
       */
      void dropFreqSeedSet()
      {
        this->mashimizerIndex.erase(
          std::remove_if(mashimizerIndex.begin(), mashimizerIndex.end(), [&] 
            (auto& mi) {return this->frequentSeeds.find(mi.hash) != this->frequentSeeds.end();}
          ), mashimizerIndex.end()
        );
      }

      public:

      int getFreqThreshold() const
      {
        return this->freqThreshold;
      }

      bool isFreqSeed(hash_t h) const
      {
        return frequentSeeds.find(h) != frequentSeeds.end();
      }


    }; //End of class Sketch
} //End of namespace skch

#endif
