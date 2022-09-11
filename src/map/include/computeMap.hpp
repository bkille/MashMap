/**
 * @file    computeMap.hpp
 * @brief   implments the sequence mapping logic
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef SKETCH_MAP_HPP 
#define SKETCH_MAP_HPP

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <zlib.h>  
#include <cassert>
#include <numeric>
#include <iostream>
//#include <cmath>

//Own includes
#include "map/include/base_types.hpp"
#include "map/include/map_parameters.hpp"
#include "map/include/commonFunc.hpp"
#include "map/include/winSketch.hpp"
#include "map/include/map_stats.hpp"
#include "map/include/slidingMap.hpp"
#include "map/include/MIIteratorL2.hpp"
#include "map/include/ThreadPool.hpp"
#include "map/include/filter.hpp"

//External includes
//#include "assert.hpp"

namespace skch
{
  /**
   * @class     skch::Map
   * @brief     L1 and L2 mapping stages
   */
  class Map
  {
    public:

      //Type for Stage L1's predicted candidate location
      struct L1_candidateLocus_t 
      {
        seqno_t seqId;                    //sequence id where read is mapped

        /* read could be mapped with its begin location
         * from [rangeStartPos, rangeEndPos]
         */
        offset_t firstPos;
        offset_t lastPos;
        int sharedSketchSize;
      };

      //Type for Stage L2's predicted mapping coordinate within each L1 candidate
      struct L2_mapLocus_t 
      {
        seqno_t seqId;                    //sequence id where read is mapped
        offset_t meanOptimalPos;          //Among multiple consecutive optimal positions, save the avg.
        //Sketch::MIIter_t optimalStart;    //optimal start mapping position (begin iterator)
        //Sketch::MIIter_t optimalEnd;      //optimal end mapping position (end iterator) 
        skch::strand_t strand;            //strand of optimal mapping
        int sharedSketchSize;             //count of shared sketch elements
      };

    private:

      //algorithm parameters
      const skch::Parameters &param;

      //reference sketch
      const skch::Sketch &refSketch;

      //Container type for saving read sketches during L1 and L2 both
      typedef Sketch::MI_Type MinVec_Type;

      typedef Sketch::MIIter_t MIIter_t;

      //Custom function for post processing the results, by default does nothing 
      typedef std::function< void(const MappingResult&) > PostProcessResultsFn_t;
      PostProcessResultsFn_t processMappingResults;

      //Container to store query sequence name and length
      //used only if one-to-one filtering is ON
      std::vector<ContigInfo> qmetadata; 

    public:

      /**
       * @brief                 constructor
       * @param[in] p           algorithm parameters
       * @param[in] refSketch   reference sketch
       * @param[in] f           optional user defined custom function to post process the reported mapping results
       */
      Map(const skch::Parameters &p, const skch::Sketch &refsketch,
          PostProcessResultsFn_t f = nullptr) :
        param(p),
        refSketch(refsketch),
        processMappingResults(f)
    {
      this->mapQuery();
    }

    private:

      /**
       * @brief   parse over sequences in query file and map each on the reference
       */
      void mapQuery()
      {
        //Count of reads mapped by us
        //Some reads are dropped because of short length
        seqno_t totalReadsPickedForMapping = 0;
        seqno_t totalReadsMapped = 0;
        seqno_t seqCounter = 0;     

        std::ofstream outstrm(param.outFileName);
        MappingResultsVector_t allReadMappings;  //Aggregate mapping results for the complete run

        //Create the thread pool 
        ThreadPool<InputSeqContainer, MapModuleOutput> threadPool( [this](InputSeqContainer* e){return mapModule(e);}, param.threads);

        for(const auto &fileName : param.querySequences)
        {
          //Open the file using kseq
          gzFile fp = gzopen(fileName.c_str(), "r");
          kseq_t *seq = kseq_init(fp);

#ifdef DEBUG
          std::cout << "INFO, skch::Map::mapQuery, mapping reads in " << fileName << std::endl;
#endif

          //size of sequence
          offset_t len;

          while ((len = kseq_read(seq)) >= 0) 
          {

            if (param.filterMode == filter::ONETOONE)
              qmetadata.push_back( ContigInfo{seq->name.s, (offset_t) seq->seq.l} );

            //Is the read too short?
            if(len < param.windowSize || len < param.kmerSize || len < param.segLength)
            {

#ifdef DEBUG
              std::cout << "WARNING, skch::Map::mapQuery, read is not long enough for mapping" << std::endl;
#endif

              seqCounter++;
              continue;
            }
            else 
            {
              totalReadsPickedForMapping++;

              //Dispatch input to thread
              threadPool.runWhenThreadAvailable(new InputSeqContainer(seq->seq.s, seq->name.s, len, seqCounter));

              //Collect output if available
              while ( threadPool.outputAvailable() )
                mapModuleHandleOutput(threadPool.popOutputWhenAvailable(), allReadMappings, totalReadsMapped, outstrm);
            }

            seqCounter++;
          } //Finish reading query input file

          //Close the input file
          kseq_destroy(seq);  
          gzclose(fp);  
        }

        //Collect remaining output objects
        while ( threadPool.running() )
          mapModuleHandleOutput(threadPool.popOutputWhenAvailable(), allReadMappings, totalReadsMapped, outstrm);

        //Filter over reference axis and report the mappings
        if (param.filterMode == filter::ONETOONE)
        {
          skch::Filter::ref::filterMappings(allReadMappings, this->refSketch);

          //Re-sort mappings by input order of query sequences
          //This order may be needed for any post analysis of output
          std::sort(allReadMappings.begin(), allReadMappings.end(), [](const MappingResult &a, const MappingResult &b)  
          {
            return (a.querySeqId < b.querySeqId);
          });

          reportReadMappings(allReadMappings, "", outstrm);
        }

        std::cout << "INFO, skch::Map::mapQuery, [count of mapped reads, reads qualified for mapping, total input reads] = [" << totalReadsMapped << ", " << totalReadsPickedForMapping << ", " << seqCounter << "]" << std::endl;
      }

      /**
       * @brief               main mapping function given an input read
       * @details             this function is run in parallel by multiple threads
       * @param[in]   input   input read details
       * @return              output object containing the mappings
       */
      MapModuleOutput* mapModule (InputSeqContainer* input)
      {
        MapModuleOutput* output = new MapModuleOutput();

        //save query sequence name
        output->qseqName = input->seqName;

        if(! param.split)   
        {
          QueryMetaData <MinVec_Type> Q;
          Q.seq = &(input->seq)[0u];
          Q.len = input->len;
          Q.seqCounter = input->seqCounter;

          MappingResultsVector_t l2Mappings;   

          //Map this sequence
          mapSingleQueryFrag(Q, l2Mappings);

          output->readMappings.insert(output->readMappings.end(), l2Mappings.begin(), l2Mappings.end());
        }
        else  //Split read mapping
        {
          int noOverlapFragmentCount = input->len / param.segLength;
          bool mappingReported = false;

          //Map individual non-overlapping fragments in the read
          for (int i = 0; i < noOverlapFragmentCount; i++)
          {
            //Prepare fragment sequence object 
            QueryMetaData <MinVec_Type> Q;
            Q.seq = &(input->seq)[0u] + i * param.segLength;
            Q.len = param.segLength;
            Q.startPos = i*param.segLength;
            Q.seqCounter = input->seqCounter;
            Q.seqName = input->seqName;

            MappingResultsVector_t l2Mappings;   

            //Map this fragment
            mapSingleQueryFrag(Q, l2Mappings);

            //Clear the seed hit table
            Q.seedHits.clear();

            //Adjust query coordinates and length in the reported mapping
            std::for_each(l2Mappings.begin(), l2Mappings.end(), [&](MappingResult &e){ 
                e.queryLen = input->len;
                e.queryStartPos = i * param.segLength;
                e.queryEndPos = i * param.segLength + Q.len - 1;
                });

            output->readMappings.insert(output->readMappings.end(), l2Mappings.begin(), l2Mappings.end());
          }

          //Map last overlapping fragment to cover the whole read
          if (noOverlapFragmentCount >= 1 && input->len % param.segLength != 0)
          {
            //Prepare fragment sequence object 
            QueryMetaData <MinVec_Type> Q;
            Q.seq = &(input->seq)[0u] + input->len - param.segLength;
            Q.len = param.segLength;
            Q.startPos = input->len - param.segLength;
            Q.seqCounter = input->seqCounter;
            Q.seqName = input->seqName;

            MappingResultsVector_t l2Mappings;   

            //Map this fragment
            mapSingleQueryFrag(Q, l2Mappings);

            //Adjust query coordinates and length in the reported mapping
            std::for_each(l2Mappings.begin(), l2Mappings.end(), [&](MappingResult &e){ 
                e.queryLen = input->len;
                e.queryStartPos = input->len - param.segLength;
                e.queryEndPos = input->len - 1;
                });

            output->readMappings.insert(output->readMappings.end(), l2Mappings.begin(), l2Mappings.end());
          }

          //merge
          mergeMappings(output->readMappings);
        }

        //filter mappings best over query sequence axis
        if(param.filterMode == filter::MAP || param.filterMode == filter::ONETOONE)
        {
          skch::Filter::query::filterMappings(output->readMappings);
        }

        //Make sure mapping boundary don't exceed sequence lengths
        this->mappingBoundarySanityCheck(input, output->readMappings);


        return output;
      }

      /**
       * @brief                       routine to handle mapModule's output of mappings
       * @param[in] output            mapping output object
       * @param[in] allReadMappings   vector to store mappings of all reads (optional use depending on filter)
       * @param[in] totalReadsMapped  counter to track count of reads mapped
       * @param[in] outstrm           outstream stream object 
       */
      template <typename Vec>
        void mapModuleHandleOutput(MapModuleOutput* output, Vec &allReadMappings, seqno_t &totalReadsMapped,
            std::ofstream &outstrm)
        {
          if(output->readMappings.size() > 0)
            totalReadsMapped++;

          if (param.filterMode == filter::ONETOONE)
          {
            //Save for another filtering round
            allReadMappings.insert(allReadMappings.end(), output->readMappings.begin(), output->readMappings.end());
          }
          else
          {  
            //Report mapping
            reportReadMappings(output->readMappings, output->qseqName, outstrm); 
          }

          delete output;
        }

      /**
       * @brief                   map the parsed query sequence (L1 and L2 mapping)
       * @param[in]   Q           metadata about query sequence
       * @param[in]   outstrm     outstream stream where mappings will be reported
       * @param[out]  l2Mappings  Mapping results in the L2 stage
       */
      template<typename Q_Info, typename VecOut>
        void mapSingleQueryFrag(Q_Info &Q, VecOut &l2Mappings)
        {
#ifdef ENABLE_TIME_PROFILE_L1_L2
          auto t0 = skch::Time::now();
#endif
          //Get seed hits
          getSeedHits(Q);

          //Mapping
          std::vector<L1_candidateLocus_t> L1_candidates;
          std::vector<skch::IntervalPoint> intervalPoints; 
          getSeedIntervalPoints(Q, intervalPoints);
          computeL1Windows(Q, intervalPoints, L1_candidates);

          //std::cout << refSketch.sortedOpenPoints.size() << ", " << refSketch.sortedClosePoints.size() << std::endl;
          //std::cout << "L1 windows...\n";
          //for (auto& l1_candidate : L1_candidates) {
            //std::cout << "Candidate window @ seqId: " << l1_candidate.seqId << ":" << l1_candidate.firstPos << ":" << l1_candidate.lastPos << std::endl;
          //}
          //std::cout << "-------------\n\n";

          std::vector<L2_mapLocus_t> temp_L2_mappings;
          //if (param.split) {
            //computeSplitReadMappedRegions(Q, refSketch.sortedOpenPoints, refSketch.sortedClosePoints, L1_candidates, temp_L2_mappings);
          //} else {
            computeFullReadMappedRegions(Q, refSketch.mashimizerIndex, L1_candidates, temp_L2_mappings);
          //}

          addL2MappingMetadata(Q, temp_L2_mappings, l2Mappings);


#ifdef ENABLE_TIME_PROFILE_L1_L2
          {
            std::chrono::duration<double> timeSpentL2 = skch::Time::now() - t1;
            std::chrono::duration<double> timeSpentMappingFragment = skch::Time::now() - t0;

            std::cerr << Q.seqCounter << " " << Q.len
              << " " << timeSpentL1.count() 
              << " " << timeSpentL2.count()
              << " " << timeSpentMappingFragment.count()
              << "\n";
          }
#endif
        }

      template <typename Q_Info>
        void getSeedHits(Q_Info &Q)
        {
          CommonFunc::sketchSequence(Q.minimizerTableQuery, Q.seq, Q.len, param.kmerSize, param.alphabetSize, param.sketchSize, Q.seqCounter);
          if(Q.minimizerTableQuery.size() == 0)
            return;

          // TODO remove them from the original sketch instead of removing for each read
          auto new_end = std::remove_if(Q.minimizerTableQuery.begin(), Q.minimizerTableQuery.end(), [&](auto mi) {
            return refSketch.isFreqSeed(mi.hash);
          });
          Q.minimizerTableQuery.erase(new_end, Q.minimizerTableQuery.end());

          Q.sketchSize = Q.minimizerTableQuery.size();
        } 

      //template <typename Q_Info, typename Vec>
        //void getL1Regions(Q_info &Q, Vec L1_regions) {


        //}

      /**
       * @brief       Find candidate regions for a read using level 1 (seed-hits) mapping
       * @details     The count of hits that should occur within a region on the reference is 
       *              determined by the threshold similarity
       *              The resulting start and end target offsets on reference is (are) an 
       *              overestimate of the mapped region. Computing better bounds is left for
       *              the following L2 stage.
       * @param[in]   Q                         query sequence details 
       * @param[out]  l1Mappings                all the read mapping locations
       */
      template <typename Q_Info, typename Vec>
        void getSeedIntervalPoints(Q_Info &Q, Vec& intervalPoints)
        {

#ifdef DEBUG
          std::cout << "INFO, skch::Map::getSeedHits, read id " << Q.seqCounter << ", minimizer count = " << Q.minimizerTableQuery.size() << " " << Q.len << "\n";
#endif

          //For invalid query (example : just NNNs), we may be left with 0 sketch size
          //Ignore the query in this case
          if(Q.minimizerTableQuery.size() == 0)
            return;

          //for (auto& mi : Q.minimizerTableQuery) 
            //std::cout << mi << std::endl;


          for(auto it = Q.minimizerTableQuery.begin(); it != Q.minimizerTableQuery.end(); it++)
          {
            //Check if hash value exists in the reference lookup index
            auto seedFind = refSketch.minimizerPosLookupIndex.find(it->hash);

            if(seedFind != refSketch.minimizerPosLookupIndex.end())
            {
              auto& hitPositionList = seedFind->second;

              // Let the strand of the hits denote wrt the reference. 
              // i.e. if - query mashimizer hits a - ref mashimizer, we mark the strand as fwd. 
              std::for_each(hitPositionList.begin(), hitPositionList.end(), [&](auto& mi) {
                  //mi.strand = mi.strand * it->strand;
                  intervalPoints.push_back(IntervalPoint {side::OPEN, mi.seqId, mi.wpos, strand_t(mi.strand * it->strand), mi.hash});
                  intervalPoints.push_back(IntervalPoint {side::CLOSE, mi.seqId, mi.wpos_end, strand_t(mi.strand * it->strand), mi.hash});
              });
            }
          }

          //Sort all the hit positions
          std::sort(intervalPoints.begin(), intervalPoints.end());

#ifdef DEBUG
          std::cout << "INFO, skch::Map:getSeedHits, read id " << Q.seqCounter << ", Count of seed hits in the reference = " << intervalPoints.size() / 2 << "\n";
#endif
          //for (auto& pi : intervalPoints) 
            //std::cout << pi << std::endl;

        }


      template <typename Q_Info, typename VecIn, typename VecOut>
        void computeL1Windows(Q_Info& Q, VecIn& intervalPoints, VecOut& L1_candidates) 
        {
#ifdef DEBUG
          std::cout << "INFO, skch::Map:computeL1Windows, read id " << Q.seqCounter << std::endl;
#endif

          int overlapCount = 0;
          int bestSketchSize = 0;

          auto trailingIt = intervalPoints.begin();
          auto leadingIt = intervalPoints.begin();

          int windowLen = Q.len - param.segLength;

          std::unordered_map<hash_t, int> hash_to_freq;

          for (auto& ip : intervalPoints) {
            overlapCount += ip.side;
            assert(overlapCount >= 0);
            assert(overlapCount <= Q.sketchSize);
          }


          while (leadingIt != intervalPoints.end())
          {
            while (
                trailingIt != intervalPoints.end() 
                && (trailingIt->seqId == leadingIt->seqId && trailingIt->pos <= leadingIt->pos - windowLen
                  || trailingIt->seqId < leadingIt->seqId))
            {
              if (trailingIt->side == side::CLOSE) {
                hash_to_freq[trailingIt->hash]--;
                if (hash_to_freq[trailingIt->hash] == 0)
                  overlapCount--;
              }
              trailingIt++;
            }
            auto currentPos = leadingIt->pos;
            while (leadingIt != intervalPoints.end() && leadingIt->pos == currentPos) {
              if (leadingIt->side == side::OPEN) {
                if (hash_to_freq[leadingIt->hash] == 0)
                  overlapCount++;
                hash_to_freq[leadingIt->hash]++;
              }
              leadingIt++;
            }
            //Is this sliding window the best we have so far?
            if (overlapCount > bestSketchSize)
            {
              bestSketchSize = overlapCount;
            }
            assert(overlapCount >= 0);
            assert(overlapCount <= Q.sketchSize);
          }

          //std::cout << "OVERLAP: " << bestSketchSize << std::endl;
          bool in_candidate = false;
          L1_candidateLocus_t l1_out = {};
          trailingIt = intervalPoints.begin();
          leadingIt = intervalPoints.begin();
          overlapCount = 0;

          while (leadingIt != intervalPoints.end())
          {
            while (
                trailingIt != intervalPoints.end() 
                && (trailingIt->seqId == leadingIt->seqId && trailingIt->pos <= leadingIt->pos - windowLen
                  || trailingIt->seqId < leadingIt->seqId))
            {
              if (trailingIt->side == side::CLOSE) {
                overlapCount--;
              }
              trailingIt++;
            }
            auto currentPos = leadingIt->pos;
            while (leadingIt != intervalPoints.end() && leadingIt->pos == currentPos) {
              if (leadingIt->side == side::OPEN) {
                overlapCount++;
              }
              leadingIt++;
            }
            IntervalPoint& ip = *leadingIt; 
            if (overlapCount > bestSketchSize * 0.9)
            {
              if (!in_candidate) {
                l1_out.firstPos = ip.pos - windowLen;
                l1_out.seqId = ip.seqId;
              }
              l1_out.lastPos = ip.pos;
              in_candidate = true;
            } else {
              if (in_candidate) {
                L1_candidates.push_back(l1_out);
              }
              in_candidate = false;
              L1_candidateLocus_t l1_out = {};
            }
          }
          if (in_candidate) {
            // Save and reset
            auto& ip = intervalPoints.back();
            l1_out.lastPos = ip.pos;
            l1_out.seqId = ip.seqId;
            L1_candidates.push_back(l1_out);
          }
#ifdef DEBUG
          std::cout << "INFO, skch::Map:computeL1Windows, Found " << L1_candidates.size() << " windows for read " << Q.seqCounter << std::endl;
#endif
        }

      /**
       * @brief                                 Revise L1 candidate regions to more precise locations
       * @param[in]   Q                         query sequence information
       * @param[in]   l1Mappings                candidate regions for query sequence found at L1
       * @param[out]  l2Mappings                Mapping results in the L2 stage
       */

      template <typename Q_Info, typename VecIn, typename VecOut>
        void addL2MappingMetadata(Q_Info &Q, VecIn &l2_vec, VecOut& l2Mappings)
        {
          for (auto& l2 : l2_vec) {
            //Compute mash distance using calculated jaccard
            //float mash_dist = Stat::j2md(1.0 * l2.sharedSketchSize/Q.sketchSize, param.kmerSize);
            float mash_dist_old = Stat::j2md(1.0 * l2.sharedSketchSize/Q.sketchSize, param.kmerSize);

            ////Correction to remove Poisson approximation from Mash Distance [BelBasi et al. (2022)]
            float mash_dist = 1 - std::exp(-mash_dist_old);

            //Compute lower bound to mash distance within 90% confidence interval
            float mash_dist_lower_bound = Stat::md_lower_bound(mash_dist, Q.sketchSize, param.kmerSize, skch::fixed::confidence_interval);

            float nucIdentity = 100 * (1 - mash_dist);
            float nucIdentityUpperBound = 100 * (1 - mash_dist_lower_bound);

            //Report the alignment
            if(nucIdentityUpperBound >= param.percentageIdentity)
            {
              MappingResult res;

              //Save the output
              {
                res.queryLen = Q.len;
                res.refStartPos = l2.meanOptimalPos ;
                res.refEndPos = l2.meanOptimalPos + Q.len - 1;
                res.queryStartPos = 0;
                res.queryEndPos = Q.len - 1;
                res.refSeqId = l2.seqId;
                res.querySeqId = Q.seqCounter;
                res.nucIdentity = nucIdentity;
                res.nucIdentityUpperBound = nucIdentityUpperBound;
                res.sketchSize = Q.sketchSize;
                res.conservedSketches = l2.sharedSketchSize;
                res.strand = l2.strand; 

                //Compute additional statistics -> strand, reference compexity
                {
                  //SlideMapper<Q_Info> slidemap(Q);
                  //slidemap.insert_ref(l2.optimalStart, l2.optimalEnd);
                  //int strandVotes, uniqueRefHashes;
                  //slidemap.computeStatistics(strandVotes, uniqueRefHashes);

                  //res.strand = strandVotes > 0 ? strnd::FWD : strnd::REV;
                  //res.strand = strnd::FWD;
                }

                l2Mappings.push_back(res);
              }
            }
          }

#ifdef DEBUG
          //std::cout << "INFO, skch::Map:doL2Mapping, read id " << Q.seqCounter << ", count of L1 candidates= " << l1Mappings.size() << ", count of L2 candidates= " << l2Mappings.size() << std::endl;
#endif
        }

      /**
       * @brief                                 Find optimal mapping within an L1 candidate
       * @param[in]   Q                         query sequence information
       * @param[in]   candidateLocus            L1 candidate location
       * @param[out]  l2_out                    L2 mapping inside L1 candidate 
       */
      template <typename Q_Info>
        void computeFullReadMappedRegions(Q_Info &Q, 
            std::vector<MashimizerInfo> mashimizerIndex,
            std::vector<L1_candidateLocus_t> &l1_vec_in,
            std::vector<L2_mapLocus_t> &l2_vec_out)
        {
//#ifdef DEBUG
          //std::cout << "INFO, skch::Map:computeL2MappedRegions, read id " << Q.seqName << "_" << Q.startPos << std::endl; 
//#endif
          for (L1_candidateLocus_t& l1_candidate : l1_vec_in) 
          {
            //std::cout << "Candidate window @ seqId: " << l1_candidate.seqId << ":" << l1_candidate.firstPos << ":" << l1_candidate.lastPos << std::endl;
           
            // Get first potential mashimizer
            const MashimizerInfo first_mashimizer = MashimizerInfo {0, l1_candidate.seqId, l1_candidate.firstPos - param.segLength, 0, 0};
            auto firstOpenIt = std::lower_bound(mashimizerIndex.begin(), mashimizerIndex.end(), first_mashimizer); 
            //std::cout << "First mashimizer @ index " << std::distance(mashimizerIndex.begin(), firstOpenIt) << std::endl;;

            // Keeps track of the lowest end position
            std::vector<skch::MashimizerInfo> slidingWindow;

            // Used to make a min-heap
            constexpr auto heap_cmp = [](const skch::MashimizerInfo& l, const skch::MashimizerInfo& r) {return l.wpos_end >= r.wpos_end;};

            // windowIt keeps track of beginning of window
            auto windowIt = firstOpenIt;

            //
            int windowLen = Q.len - param.segLength;
            
            // slideMap tracks the S(A or B) and S(A) and S(B)
            SlideMapper<Q_Info> slideMap(Q);

            int beginOptimalPos = 0;
            int lastOptimalPos = 0;
            int bestSketchSize = 0;
            bool in_candidate = false;
            L2_mapLocus_t l2_out = {};

            // Set up the window
            while (windowIt != mashimizerIndex.end() && windowIt->seqId == l1_candidate.seqId && windowIt->wpos < l1_candidate.firstPos) 
            {
              if (windowIt->wpos_end > l1_candidate.firstPos) {
                  slidingWindow.push_back(*windowIt);
                  std::push_heap(slidingWindow.begin(), slidingWindow.end(), heap_cmp);
                  slideMap.insert_mashimizer(*windowIt);
              }
              windowIt++;
            }
            //std::cout << "Start mashimizers loaded\n";

            while (windowIt->seqId == l1_candidate.seqId && windowIt->wpos < l1_candidate.lastPos) 
            {
              int prev_strand_votes = slideMap.strand_votes;
              while (!slidingWindow.empty() && slidingWindow.front().wpos_end <= windowIt->wpos - windowLen) {
                // Remove mashimizer from  sorted window
                slideMap.delete_mashimizer(slidingWindow.front());

                // Remove mashimizer from end-ordered heap
                std::pop_heap(slidingWindow.begin(), slidingWindow.end(), heap_cmp);
                slidingWindow.pop_back();
              }
              slideMap.insert_mashimizer(*windowIt);
              slidingWindow.push_back(*windowIt);
              std::push_heap(slidingWindow.begin(), slidingWindow.end(), heap_cmp);

              //std::cout << "Leading mashimizer @ " << windowIt->seqId << ":" <<  windowIt->wpos << std::endl;
              //std::cout << "Overlap: " << slideMap.sharedSketchElements << std::endl;

              //Is this sliding window the best we have so far?
              if (slideMap.sharedSketchElements > bestSketchSize)
              {
                in_candidate = true;
                bestSketchSize = slideMap.sharedSketchElements;
                l2_out.sharedSketchSize = slideMap.sharedSketchElements;

                //Save the position
                beginOptimalPos = windowIt->wpos - windowLen;

                // TODO This should be halfway between windowIt and next(windowIt)
                lastOptimalPos = windowIt->wpos - windowLen;
              }
              else if(slideMap.sharedSketchElements == bestSketchSize)
              {
                if (!in_candidate) {
                  l2_out.sharedSketchSize = slideMap.sharedSketchElements;

                  //Save the position
                  beginOptimalPos = windowIt->wpos - windowLen;
                }

                in_candidate = true;
                //Still save the position
                lastOptimalPos = windowIt->wpos - windowLen;
              } else {
                if (in_candidate) {
                  // Save and reset
                  lastOptimalPos = windowIt->wpos - windowLen;
                  l2_out.meanOptimalPos =  (beginOptimalPos + lastOptimalPos) / 2;
                  l2_out.seqId = windowIt->seqId;
                  l2_out.strand = prev_strand_votes >= 0 ? strnd::FWD : strnd::REV;
                  l2_vec_out.push_back(l2_out);
                  l2_out = L2_mapLocus_t();
                }
                in_candidate = false;
              }
              windowIt++;
            }
            if (in_candidate) {
              // Save and reset
              //lastOptimalPos = std::prev(windowIt->wpos;
              l2_out.meanOptimalPos =  (beginOptimalPos + lastOptimalPos) / 2;
              l2_out.seqId = std::prev(windowIt)->seqId;
              l2_out.strand = slideMap.strand_votes >= 0 ? strnd::FWD : strnd::REV;
              l2_vec_out.push_back(l2_out);
            }
          }
        }

      /**
       * @brief                                 Find optimal mapping within an L1 candidate
       * @param[in]   Q                         query sequence information
       * @param[in]   candidateLocus            L1 candidate location
       * @param[out]  l2_out                    L2 mapping inside L1 candidate 
       */
      template <typename Q_Info>
        void computeSplitReadMappedRegions(Q_Info &Q, 
            std::vector<skch::IntervalPoint> &intervalPoints,
            std::vector<L2_mapLocus_t> &l2_vec_out)
        {
//#ifdef DEBUG
          //std::cout << "INFO, skch::Map:computeL2MappedRegions, read id " << Q.seqName << "_" << Q.startPos << std::endl; 
//#endif
          int strand_votes = 0;
          int overlapCount = 0;
          int beginOptimalPos = 0;
          int lastOptimalPos = 0;
          int bestSketchSize = 0;
          bool in_candidate = false;
          L2_mapLocus_t l2_out = {};

          for (int i = 0; i < intervalPoints.size(); i++) 
          {
            const IntervalPoint& ip = intervalPoints[i];
            overlapCount += ip.side;
            
            if (i == intervalPoints.size() - 1 || ip.pos != intervalPoints[i+1].pos) 
            {
              //Is this sliding window the best we have so far?
              if (overlapCount > bestSketchSize)
              {
                in_candidate = true;
                bestSketchSize = overlapCount;
                l2_out.sharedSketchSize = overlapCount;

                //Save the position
                beginOptimalPos = ip.pos;
                lastOptimalPos = ip.pos;
              }
              else if(overlapCount == bestSketchSize)
              {
                if (!in_candidate) {
                  l2_out.sharedSketchSize = overlapCount;

                  //Save the position
                  beginOptimalPos = ip.pos;
                }

                in_candidate = true;
                //Still save the position
                lastOptimalPos = ip.pos;
              } else {
                if (in_candidate) {
                  // Save and reset
                  lastOptimalPos = ip.pos;
                  l2_out.meanOptimalPos =  (beginOptimalPos + lastOptimalPos) / 2;
                  l2_out.seqId = ip.seqId;
                  l2_out.strand = strand_votes >= 0 ? strnd::FWD : strnd::REV;
                  l2_vec_out.push_back(l2_out);
                  l2_out = L2_mapLocus_t();
                }
                in_candidate = false;
              }
            }

            strand_votes += ip.strand * ip.side;
            //std::cout << overlapCount << ", " << strand_votes << "\t" << (ip.side == side::OPEN ? "OPEN " : "CLOSE") << " @ " << std::to_string(ip.strand) << ", " << ip.seqId << ":" << ip.pos << "\t" << std::endl;
            assert(overlapCount >= 0);
            assert(overlapCount <= Q.sketchSize);
          }//End of while loop

          if (in_candidate) {
            // Save and reset
            auto& ip = intervalPoints.back();
            lastOptimalPos = ip.pos;
            l2_out.meanOptimalPos =  (beginOptimalPos + lastOptimalPos) / 2;
            l2_out.seqId = ip.seqId;
            l2_out.strand = strand_votes >= 0 ? strnd::FWD : strnd::REV;
            l2_vec_out.push_back(l2_out);
          }

        }

        struct compareMinimizersByPos
        {
          typedef std::pair<seqno_t, offset_t> P;

          bool operator() (const MinimizerInfo &m, const P &val)
          {
            return ( P(m.seqId, m.wpos) < val);
          }

          bool operator() (const P &val, const MinimizerInfo &m)
          {
            return (val < P(m.seqId, m.wpos) );
          }
        } cmp;

        MIIter_t searchIndex(Sketch::MI_Type& mindex, seqno_t seqId, offset_t winpos) const
        {
            std::pair<seqno_t, offset_t> searchPosInfo(seqId, winpos);

            /*
            *          * std::lower_bound --  Returns an iterator pointing to the first element in the range
            *                   *                      that is not less than (i.e. greater or equal to) value.
            *                            */
            MIIter_t iter = std::lower_bound(mindex.begin(), mindex.end(), searchPosInfo, this->cmp);

            return iter;
        }////Define std::map and let it contain only the query minimizers

      /**
       * @brief                                 Compute jaccard similarity of given query seq 
       *                                        at a single position on reference sequence
       * @param[in]   Q                         query sequence information
       * @param[in]   seqId                     reference sequence id
       * @param[in]   refStartPos               offset on reference sequence
       * return                                 jaccard similarity
       */
      //template <typename Q_Info>
        //double computeJaccardSinglePos(Q_Info &Q, seqno_t seqId, offset_t refStartPos)
        //{
          ////Look up L1 candidate's begin in the index
          //MIIter_t superWindowRangeStart = this->refSketch.searchIndex(seqId, refStartPos);

          //if (Q.sketchSize == 0)
            //return 0;
          
          ////If iterator points to index end or a different reference sequence, return zero 
          //if ( this->refSketch.isMinimizerIndexEnd(superWindowRangeStart) || superWindowRangeStart->seqId != seqId)
            //return 0;

          ////Count of minimizer windows in a super-window
          //offset_t countMinimizerWindows = Q.len - (param.windowSize-1) - (param.kmerSize-1); 

          ////Look up the end of the first L2 super-window in the index
          //MIIter_t superWindowRangeEnd = this->refSketch.searchIndex(seqId, 
              //superWindowRangeStart->wpos + countMinimizerWindows);

          //SlideMapper<Q_Info> slidemap(Q);

          ////Insert all the minimizers in the first super-window
          //slidemap.insert_ref(superWindowRangeStart, superWindowRangeEnd);

          //return slidemap.sharedSketchElements * 1.0 / Q.sketchSize;
        //}


      /**
       * @brief                       Merge the consecutive fragment mappings reported in each query 
       * @param[in/out] readMappings  Mappings computed by Mashmap (L2 stage) for a read
       */
      template <typename VecIn>
        void mergeMappings(VecIn &readMappings)
        {
          assert(param.split == true);

          if(readMappings.size() < 2)
            return;

          // Split reads must be nearly adjacent in order to be chained
          uint32_t cutoffDistance = 3.0 * (param.segLength / (double)param.sketchSize);

          //Sort the mappings by reference position
          std::sort(readMappings.begin(), readMappings.end(), [](const MappingResult &a, const MappingResult &b)  
              {
              return std::tie(a.refSeqId, a.refStartPos, a.queryStartPos) < std::tie(b.refSeqId, b.refStartPos, b.queryStartPos);
              });

          //First assign a unique id to each split mapping in the sorted order
          for(auto it = readMappings.begin(); it != readMappings.end(); it++)
          {
            it->splitMappingId = std::distance(readMappings.begin(), it);
          }

          // Keep track of best sketches
          std::vector<float> bestSketchSizeProp(std::ceil(readMappings[0].queryLen / (float)param.segLength));

          //Start the procedure to identify the chains
          for(auto it = readMappings.begin(); it != readMappings.end(); it++)
          {
            //Which fragment is this wrt. the complete read
            uint32_t currMappingFragno = std::ceil(it->queryStartPos * 1.0/param.segLength);

            if (static_cast<float>(it->conservedSketches) / it->sketchSize > bestSketchSizeProp[currMappingFragno]) {
              bestSketchSizeProp[currMappingFragno] = static_cast<float>(it->conservedSketches) / it->sketchSize;
            }
          }
          for(auto it = readMappings.begin(); it != readMappings.end(); it++)
          {
            uint32_t currMappingFragno = std::ceil(it->queryStartPos * 1.0/param.segLength);
            if (static_cast<float>(it->conservedSketches) / it->sketchSize < .99*bestSketchSizeProp[currMappingFragno]) {
              //it->discard = 1;
              continue;
            }
            for(auto it2 = std::next(it); it2 != readMappings.end(); it2++)
            {
              auto thisMappingFragno = std::ceil(it2->queryStartPos * 1.0/ param.segLength);
              if (static_cast<float>(it2->conservedSketches) / it->sketchSize  < .99*bestSketchSizeProp[thisMappingFragno]) {
                //it->discard = 1;
                continue;
              }

              // If mappings are on different seqs or too far apart, stop looking
              if(it2->refSeqId != it->refSeqId || std::abs(it2->refStartPos - it->refEndPos) > param.segLength + 1)
                  break;
             
              // Check if it is consecutive query fragment and strand matches
              if( it2->strand == it->strand
                  && thisMappingFragno == currMappingFragno + (it->strand == strnd::FWD ? 1 : -1) )
              {

                // Use offset to identify appropriate chaining for end-of-read fragments.
                int offset = 0; //it->queryEndPos - it2->queryStartPos + 1;
                if (it->strand == strnd::FWD) {
                  offset = it->queryEndPos - it2->queryStartPos + 1;
                } else {
                  offset = it2->queryEndPos - it->queryStartPos + 1;
                }

                //If this mapping is too far from current mapping being evaluated, skip
                if (std::abs(it2->refStartPos + offset - it->refEndPos) > cutoffDistance)
                  continue;

                it2->splitMappingId = it->splitMappingId;   //merge
              }
            }
          }

          //Sort the mappings by post-merge split mapping id
          std::sort(readMappings.begin(), readMappings.end(), [](const MappingResult &a, const MappingResult &b)  
              {
              return a.splitMappingId < b.splitMappingId;
              });

          if (!param.noChain) {
            for(auto it = readMappings.begin(); it != readMappings.end();)
            {
              //Bucket by each chain
              auto it_end = std::find_if(it, readMappings.end(), [&](const MappingResult &e){return e.splitMappingId != it->splitMappingId;} ); 

              //[it -- it_end) represents same chain


              //Incorporate chain information into first mapping

              //compute chain length
              std::for_each(it, it_end, [&](MappingResult &e)
              {
                // Only use the best mappings in a chain to compute the window
                it->queryStartPos = std::min( it->queryStartPos, e.queryStartPos);
                it->refStartPos = std::min( it->refStartPos, e.refStartPos);

                it->queryEndPos = std::max( it->queryEndPos, e.queryEndPos);
                it->refEndPos = std::max( it->refEndPos, e.refEndPos);
              });

              //Mean identity of all mappings in the chain
              //TODO some of these overlap, leading to inaccurate averages
              it->nucIdentity = (   std::accumulate(it, it_end, 0.0, 
                                    [](double x, MappingResult &e){ return x + e.nucIdentity; })     
                  )/ std::distance(it, it_end);

              //Discard other mappings of this chain
              std::for_each( std::next(it), it_end, [&](MappingResult &e){ e.discard = 1; });

              //advance the iterator
              it = it_end;
            }

            //for (auto& m : readMappings) {
              //std::cout << "Final mapping with " << m.conservedSketches << " shared sketcheds" << std::endl;
            //}
            readMappings.erase( 
                std::remove_if(readMappings.begin(), readMappings.end(), [&](MappingResult &e){ return e.discard == 1; }),
                readMappings.end());
            //std::cout << "AFTER\n";
            //for (auto& m : readMappings) {
              //std::cout << "Final mapping with " << m.conservedSketches << " shared sketcheds" << std::endl;
            //}
          }
       }

      /**
       * @brief                       This routine is to make sure that all mapping boundaries
       *                              on query and reference are not outside total 
       *                              length of sequeunces involved
       * @param[in]     input         input read details
       * @param[in/out] readMappings  Mappings computed by Mashmap (L2 stage) for a read
       */
      template <typename VecIn>
        void mappingBoundarySanityCheck(InputSeqContainer* input, VecIn &readMappings)
        {
          for(auto &e : readMappings)
          {
            //reference start pos
            {
              if(e.refStartPos < 0)
                e.refStartPos = 0;
              if(e.refStartPos >= this->refSketch.metadata[e.refSeqId].len)
                e.refStartPos = this->refSketch.metadata[e.refSeqId].len - 1;
            }

            //reference end pos
            {
              if(e.refEndPos < e.refStartPos)
                e.refEndPos = e.refStartPos;
              if(e.refEndPos >= this->refSketch.metadata[e.refSeqId].len)
                e.refEndPos = this->refSketch.metadata[e.refSeqId].len - 1;
            }

            //query start pos
            {
              if(e.queryStartPos < 0)
                e.queryStartPos = 0;
              if(e.queryStartPos >= input->len)
                e.queryStartPos = input->len;
            }

            //query end pos
            {
              if(e.queryEndPos < e.queryStartPos)
                e.queryEndPos = e.queryStartPos;
              if(e.queryEndPos >= input->len)
                e.queryEndPos = input->len;
            }
          }
        }

      /**
       * @brief                         Report the final read mappings to output stream
       * @param[in]   readMappings      mapping results for single or multiple reads
       * @param[in]   queryName         input required if reporting one read at a time
       * @param[in]   outstrm           file output stream object
       */
      void reportReadMappings(MappingResultsVector_t &readMappings, const std::string &queryName, 
          std::ofstream &outstrm)
      {
        //Print the results
        for(auto &e : readMappings)
        {
          assert(e.refSeqId < this->refSketch.metadata.size());

          outstrm  << (param.filterMode == filter::ONETOONE ? qmetadata[e.querySeqId].name : queryName)
            << " " << e.queryLen 
            << " " << e.queryStartPos
            << " " << e.queryEndPos
            << " " << (e.strand == strnd::FWD ? "+" : "-") 
            << " " << this->refSketch.metadata[e.refSeqId].name
            << " " << this->refSketch.metadata[e.refSeqId].len
            << " " << e.refStartPos 
            << " " << e.refEndPos
            << " " << e.nucIdentity
            << " " << e.conservedSketches;

#ifdef DEBUG
          outstrm << std::endl;
#else
          outstrm << "\n";
#endif

          //User defined processing of the results
          if(processMappingResults != nullptr)
            processMappingResults(e);
        }
      }

    public:

      /**
       * @brief     An optional utility function to save the 
       *            reported results by the L2 stage into a vector
       */
      static void insertL2ResultsToVec(MappingResultsVector_t &v, const MappingResult &reportedL2Result)
      {
        v.push_back(reportedL2Result);
      }

  };

}

#endif
