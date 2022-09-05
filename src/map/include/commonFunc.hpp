/**
 * @file    commonFunc.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef COMMON_FUNC_HPP 
#define COMMON_FUNC_HPP

#include <unordered_set>
#include <utility>
#include <vector>
#include <algorithm>
#include <deque>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <cassert>
#include <memory>

//Own includes
#include "map/include/map_parameters.hpp"
#include "map/include/base_types.hpp"

//External includes
#include "common/murmur3.h"
#include "common/prettyprint.hpp"
#include "assert.hpp"

namespace skch
{
  /**
   * @namespace skch::CommonFunc
   * @brief     Implements frequently used common functions
   */
  namespace CommonFunc
  {

    template <typename I>
    struct Pivot {
        I p;
        int64_t rank;
    };

    //seed for murmerhash
    const int seed = 42;

    /**
     * @brief   reverse complement of kmer (borrowed from mash)
     * @note    assumes dest is pre-allocated
     */
    inline void reverseComplement(const char* src, char* dest, int length) 
    {
      for ( int i = 0; i < length; i++ )
      {    
        char base = src[i];

        switch ( base )
        {    
          case 'A': base = 'T'; break;
          case 'C': base = 'G'; break;
          case 'G': base = 'C'; break;
          case 'T': base = 'A'; break;
          default: break;
        }    

        dest[length - i - 1] = base;
      }    
    }

    /**
     * @brief               convert DNA or AA alphabets to upper case
     * @param[in]   seq     pointer to input sequence
     * @param[in]   len     length of input sequence
     */
    inline void makeUpperCase(char* seq, offset_t len)
    {
      for ( int i = 0; i < len; i++ )
      {
        if (seq[i] > 96 && seq[i] < 123)
        {
          seq[i] -= 32;
        }
      }
    }

    /**
     * @brief   hashing kmer string (borrowed from mash)
     */
    inline hash_t getHash(const char* seq, int length)
    {
      char data[16];
      MurmurHash3_x64_128(seq, length, seed, data);

      hash_t hash;

      hash = *((hash_t *)data);

      return hash;
    }

    /**
     * @brief       Compute the minimum s kmers for a string.
     * @param[out]  mashimizerIndex     container storing sketched Kmers 
     * @param[in]   seq                 pointer to input sequence
     * @param[in]   len                 length of input sequence
     * @param[in]   kmerSize
     * @param[in]   s                   sketch size. 
     * @param[in]   seqCounter          current sequence number, used while saving the position of minimizer
     */
    template <typename T>
      inline void sketchSequence(
          std::vector<T> &mashimizerIndex, 
          char* seq, 
          offset_t len,
          int kmerSize, 
          int alphabetSize,
          int sketchSize,
          seqno_t seqCounter)
    {
        makeUpperCase(seq, len);

        //Compute reverse complement of seq
        std::unique_ptr<char[]> seqRev(new char[len]);
        //char* seqRev = new char[len];

        if(alphabetSize == 4) //not protein
          CommonFunc::reverseComplement(seq, seqRev.get(), len);

        // TODO cleanup
        std::unordered_map<hash_t, MashimizerInfo> sketched_vals;
        std::vector<hash_t> sketched_heap;

        for(offset_t i = 0; i < len - kmerSize + 1; i++)
        {
          //Hash kmers
          hash_t hashFwd = CommonFunc::getHash(seq + i, kmerSize); 
          hash_t hashBwd;

          if(alphabetSize == 4)
            hashBwd = CommonFunc::getHash(seqRev.get() + len - i - kmerSize, kmerSize);
          else  //proteins
            hashBwd = std::numeric_limits<hash_t>::max();   //Pick a dummy high value so that it is ignored later

          //Consider non-symmetric kmers only
          if(hashBwd != hashFwd)
          {
            //Take minimum value of kmer and its reverse complement
            hash_t currentKmer = std::min(hashFwd, hashBwd);

            //Check the strand of this minimizer hash value
            auto currentStrand = hashFwd < hashBwd ? strnd::FWD : strnd::REV;

            ////std::cout << seqCounter << "\t" << std::string(seq + i, kmerSize) << "-->" <<  currentKmer << std::endl;
            if (sketched_vals.empty() || sketched_vals.find(currentKmer) == sketched_vals.end()) {

              // Add current hash to heap
              if (sketched_vals.size() < sketchSize || currentKmer < sketched_heap[0])  {
                  sketched_vals[currentKmer] = MashimizerInfo{currentKmer, seqCounter, i, i, currentStrand};
                  sketched_heap.push_back(currentKmer);
                  std::push_heap(sketched_heap.begin(), sketched_heap.end());
              }

              // Remove one if too large
              if (sketched_vals.size() > sketchSize) {
                  sketched_vals.erase(sketched_heap[0]);
                  std::pop_heap(sketched_heap.begin(), sketched_heap.end());
                  sketched_heap.pop_back();
              }
            } else {
              // TODO these sketched values might never be useful, might save memory by deleting
              // extend the length of the window
              sketched_vals[currentKmer].wpos_end = i;
              sketched_vals[currentKmer].strand += currentStrand == strnd::FWD ? 1 : -1;
            }
          }
        }
        DEBUG_ASSERT(sketched_vals.size() <= sketchSize);
        mashimizerIndex.reserve(sketched_vals.size());
        std::for_each(sketched_vals.begin(), sketched_vals.end(),
            [&mashimizerIndex](auto& pair) {
            pair.second.strand = pair.second.strand > 0 ? strnd::FWD : (pair.second.strand == 0 ? strnd::AMBIG : strnd::REV);
            DEBUG_ASSERT(std::abs(pair.second.strand) < 2);
            mashimizerIndex.push_back(pair.second);
        });

        return;
    }
    
    /**
     * @brief       Compute winnowed mashimizers from a given sequence and add to the index
     * @param[out]  mashimizerIndex  table storing mashimizers and their position as we compute them
     * @param[in]   seq             pointer to input sequence
     * @param[in]   len             length of input sequence
     * @param[in]   kmerSize
     * @param[in]   windowSize
     * @param[in]   s               sketch size. 
     * @param[in]   seqCounter      current sequence number, used while saving the position of minimizer
     */
    template <typename T>
      inline void addMashimizers(std::vector<T> &mashimizerIndex, 
          char* seq, offset_t len,
          int kmerSize, 
          int windowSize,
          int alphabetSize,
          int sketchSize,
          seqno_t seqCounter)
      {
        /**
         * Double-ended queue (saves minimum at front end)
         * Saves pair of the minimizer and the position of hashed kmer in the sequence
         * Position of kmer is required to discard kmers that fall out of current window
         */


        std::deque< std::tuple<hash_t, strand_t, offset_t> > Q;
        using windowMap_t = std::map<hash_t, std::pair<MashimizerInfo, uint64_t>>;
        windowMap_t sortedWindow;
        Pivot<typename windowMap_t::iterator>  piv = {sortedWindow.begin(), 0};

        makeUpperCase(seq, len);

        //Compute reverse complement of seq
        std::unique_ptr<char[]> seqRev(new char[len]);

        if(alphabetSize == 4) //not protein
          CommonFunc::reverseComplement(seq, seqRev.get(), len);

        for(offset_t i = 0; i < len - kmerSize + 1; i++)
        {
          //The serial number of current sliding window
          //First valid window appears when i = windowSize - 1
          offset_t currentWindowId = i + kmerSize - windowSize;

          if (currentWindowId == 0) {
            uint64_t rank = 1;
            auto iter = sortedWindow.begin();
            while (iter != sortedWindow.end() && rank <= sketchSize) {
              iter->second.first.wpos = currentWindowId; 
              std::advance(iter, 1);
              rank += 1;
            }
          }

          //Hash kmers
          hash_t hashFwd = CommonFunc::getHash(seq + i, kmerSize); 
          hash_t hashBwd;

          if(alphabetSize == 4)
            hashBwd = CommonFunc::getHash(seqRev.get() + len - i - kmerSize, kmerSize);
          else  //proteins
            hashBwd = std::numeric_limits<hash_t>::max();   //Pick a dummy high value so that it is ignored later

          //Take minimum value of kmer and its reverse complement
          hash_t currentKmer = std::min(hashFwd, hashBwd);

          //Check the strand of this minimizer hash value
          auto currentStrand = hashFwd < hashBwd ? strnd::FWD : strnd::REV;


          //std::cout << std::endl << currentWindowId << "\t" << currentKmer << std::endl;
          //std::cout << i << "\t" << std::string(seq + i, kmerSize) << "-->" <<  currentKmer << std::endl;
          //std::cout << "Window size:\t" << sortedWindow.size() << "\tPiv:\t" << (piv.p == sortedWindow.end() ? "INFTY" : std::to_string(piv.p->first)) << std::endl;
          //for (auto it = sortedWindow.begin(); it != piv.p; it++) {
            //std::cout << "(" << it->first << ": " << it->second.second << "), ";
          //}
          //std::cout << std::endl;
          DEBUG_ASSERT(std::distance(sortedWindow.begin(), piv.p) == std::min<int>(sortedWindow.size(), sketchSize), seqCounter, currentWindowId, i);

          //If front minimum is not in the current window, remove it
          if (!Q.empty() && std::get<2>(Q.front()) <  currentWindowId) {
            const auto [leaving_hash, leaving_strand, _] = Q.front();

            // Check if we've deleted the hash already
            if (sortedWindow.find(leaving_hash) != sortedWindow.end()) {
              // If the hash that is getting popped off is still in the window and it is now leaving the window 
              // wpos != -1 and wpos_end == -1 --> still in window
              if (sortedWindow[leaving_hash].first.wpos != -1 and sortedWindow[leaving_hash].first.wpos_end == -1 && sortedWindow[leaving_hash].second == 1) {
                sortedWindow[leaving_hash].first.wpos_end = currentWindowId;
                //sortedWindow[leaving_hash].first.strand = sortedWindow[leaving_hash].first.strand >= 0 ? strnd::FWD : strnd::REV;
                //std::cout << leaving_hash << " @ rank " << std::distance(sortedWindow.begin(), sortedWindow.find(leaving_hash)) << " is leaving the window and being added to the index\n";
                mashimizerIndex.push_back(sortedWindow[leaving_hash].first);
              } else if (sortedWindow[leaving_hash].second != 1) {
                //std::cout << leaving_hash << " @ rank " << std::distance(sortedWindow.begin(), sortedWindow.find(leaving_hash)) << " is leaving the window but there are more\n";
              }

              // Remove hash
              sortedWindow[leaving_hash].second -= 1;
              if (sortedWindow[leaving_hash].second == 0) {
                if (leaving_hash == piv.p->first) {
                  std::advance(piv.p, 1);
                }
                else if (leaving_hash < piv.p->first) {
                  // Kicking out a sketched element
                  if (sortedWindow.size() >= sketchSize + 1) {
                      std::advance(piv.p, 1);
                  }
                }
                sortedWindow.erase(leaving_hash);
              } else {
                // Not removing hash, but need to adjust the strand
                if ((sortedWindow[leaving_hash].first.strand == 0 || sortedWindow[leaving_hash].first.strand - leaving_strand == 0)
                    && leaving_hash < piv.p->first) {
                  sortedWindow[leaving_hash].first.wpos_end = currentWindowId;
                  mashimizerIndex.push_back(sortedWindow[leaving_hash].first);
                  sortedWindow[leaving_hash].first.wpos = currentWindowId;
                  sortedWindow[leaving_hash].first.wpos_end = -1;
                }
                sortedWindow[leaving_hash].first.strand -= leaving_strand;
              }
            }
            Q.pop_front();
          }

          DEBUG_ASSERT(std::distance(sortedWindow.begin(), piv.p) == std::min<int>(sortedWindow.size(), sketchSize), seqCounter, currentWindowId, i);
          //Consider non-symmetric kmers only
          if(hashBwd != hashFwd)
          {
            // Add current hash to window
            Q.push_back(std::make_tuple(currentKmer, currentStrand, i)); 
            if (sortedWindow[currentKmer].second == 0) {
                auto mi = MashimizerInfo{currentKmer, seqCounter, -1, -1, currentStrand};
                sortedWindow[currentKmer].first = mi;

                if (sortedWindow.size() >= sketchSize + 2 && currentKmer < piv.p->first) {
                    piv.p--;
                }
                if (sortedWindow.size() <= sketchSize) {
                    piv.p = sortedWindow.end();
                } else if (sortedWindow.size() == sketchSize + 1) {
                    piv.p = std::prev(sortedWindow.end());
                }
            } else {
              if ((sortedWindow[currentKmer].first.strand + currentStrand == 0 || sortedWindow[currentKmer].first.strand == 0) 
                  && currentKmer < piv.p->first) {
                sortedWindow[currentKmer].first.wpos_end = currentWindowId;
                mashimizerIndex.push_back(sortedWindow[currentKmer].first);
                sortedWindow[currentKmer].first.wpos = currentWindowId;
                sortedWindow[currentKmer].first.wpos_end = -1;
              }
              sortedWindow[currentKmer].first.strand += currentStrand;
            }
            sortedWindow[currentKmer].second += 1;
          }

          DEBUG_ASSERT(std::distance(sortedWindow.begin(), piv.p) == std::min<int>(sortedWindow.size(), sketchSize), seqCounter, currentWindowId, i);

          //Select the minimizer from Q and put into index
          if(currentWindowId >= 0)
          {

            // Does the new kmer belong in the sketch?
            if (
                hashBwd != hashFwd                                                  // Non-symmetric 
                && ((piv.p == sortedWindow.end()) || (currentKmer < piv.p->first))  // Belongs in sketch
                && sortedWindow[currentKmer].first.wpos == -1) {                    // Haven't seen it in the window yet
              //std::cout << "Adding currentKmer = " << currentKmer << " to the sketch\n";
              sortedWindow[currentKmer].first.wpos = currentWindowId;
            }

            // Did we incorporate a previously hashed kmer into the sketch?
            auto& sth_mi = std::prev(piv.p)->second.first;
            if (sth_mi.wpos == -1) {
              //std::cout << "Adding bordered kmer = " << sth_mi.hash << " to the sketch\n";
              sth_mi.wpos = currentWindowId;
            }

            // Did we kick a mashimizer into non-sketch territory?
            if (piv.p != sortedWindow.end()) {
              auto& splus1th_mi = piv.p->second.first;
              if (splus1th_mi.wpos != -1) {
                //std::cout << "Removing bordered kmer = " << splus1th_mi.hash << " from the sketch and adding to index\n";
                splus1th_mi.wpos_end = currentWindowId;
                mashimizerIndex.push_back(
                    MashimizerInfo(splus1th_mi)
                );
                splus1th_mi.wpos = -1;
                splus1th_mi.wpos_end = -1;
              }
            }
#ifdef DEBUG
            DEBUG_ASSERT(std::distance(sortedWindow.begin(), piv.p) == std::min<int>(sortedWindow.size(), sketchSize), seqCounter, currentWindowId, i);
            DEBUG_ASSERT(piv.p == sortedWindow.end() || (piv.p->second.first.wpos == -1 && piv.p->second.first.wpos_end == -1));
            DEBUG_ASSERT((sortedWindow.size() == 0 || currentWindowId < 0) || (std::prev(piv.p)->second.first.wpos != -1 && std::prev(piv.p)->second.first.wpos_end == -1));
            for (auto it = sortedWindow.begin(); it != sortedWindow.end(); it++) {
              if (piv.p == sortedWindow.end() || it->first < piv.p->first) {
                DEBUG_ASSERT(it->second.first.wpos != -1, it->second.first);
                DEBUG_ASSERT(it->second.first.wpos_end == -1);
              } else {
                DEBUG_ASSERT(it->second.first.wpos == -1, it->second.first, currentWindowId);
                DEBUG_ASSERT(it->second.first.wpos_end == -1);
              }
            }
#endif
          } else {
              if (hashBwd != hashFwd && sortedWindow[currentKmer].second == 1) {
                  // Seeing kmer for the first time
                  if (sortedWindow.size() < sketchSize + 1) {
                    piv.p = sortedWindow.end();
                  } else if (sortedWindow.size() == sketchSize + 1) {
                    piv.p = std::prev(sortedWindow.end());
                  }
              }
          }

          //if (sortedWindow.size() > sketchSize*2 + 20) {
            //sortedWindow.erase(std::prev(sortedWindow.end()));
          //}
          DEBUG_ASSERT(std::distance(sortedWindow.begin(), piv.p) == std::min<int>(sortedWindow.size(), sketchSize), seqCounter, currentWindowId, i);
          DEBUG_ASSERT(piv.p == sortedWindow.end() || (piv.p->second.first.wpos == -1 && piv.p->second.first.wpos_end == -1));
          DEBUG_ASSERT((sortedWindow.size() == 0 || currentWindowId < 0) || (std::prev(piv.p)->second.first.wpos != -1 && std::prev(piv.p)->second.first.wpos_end == -1));
        }

        uint64_t rank = 1;
        auto iter = sortedWindow.begin();
        while (iter != sortedWindow.end() && rank <= sketchSize) {
          if (iter->second.first.wpos != -1) {
            iter->second.first.wpos_end = len - kmerSize + 1;
            mashimizerIndex.push_back(iter->second.first);
          }
          std::advance(iter, 1);
          rank += 1;
        }
        std::for_each(mashimizerIndex.begin(), mashimizerIndex.end(), [] (auto& mi) {
          mi.strand = mi.strand < 0 ? (mi.strand == 0 ? strnd::AMBIG : strnd::REV) : strnd::FWD;
        });
        std::sort(mashimizerIndex.begin(), mashimizerIndex.end(), [](auto& l, auto& r) {return l.wpos < r.wpos;});

        // TODO assert that this is false
        mashimizerIndex.erase(std::unique(mashimizerIndex.begin(), mashimizerIndex.end(), 
                    [](auto& l, auto& r) {
                        return (l.wpos == r.wpos) && (l.hash == r.hash);
            }),
            mashimizerIndex.end());
#ifdef DEBUG
        ////std::cout << "INFO, skch::CommonFunc::addMinimizers, inserted minimizers for sequence id = " << seqCounter << "\n";
        ////std::cout << "INFO, skch::CommonFunc::addMinimizers, length of sequence  = " << len << "\n";
        //assert(std::all_of(mashimizerIndex.begin(), mashimizerIndex.end(), [](auto& mi) {return mi.wpos >= 0;}));
        //assert(std::all_of(mashimizerIndex.begin(), mashimizerIndex.end(), [](auto& mi) {return mi.wpos_end >= 0;}));
        //std::vector<MashimizerInfo> endpos_heap;
        //auto heap_cmp = [](auto& l, auto& r) {return l.wpos_end >= r.wpos_end;};
        //for (auto& mi : mashimizerIndex) {
          //while (!endpos_heap.empty() && endpos_heap.front().wpos_end <= mi.wpos) {
            //std::pop_heap(endpos_heap.begin(), endpos_heap.end(), heap_cmp); 
            //endpos_heap.pop_back();
          //}
          //endpos_heap.push_back(mi);
          //std::push_heap(endpos_heap.begin(), endpos_heap.end(), heap_cmp);
          //assert(endpos_heap.size() <= sketchSize);
        //}
#endif
        ////std::cout << "BEFORE " << mashimizerIndex.size() << "\n";
        ////std::cout << "AFTER " << mashimizerIndex.size() << "\n";
      }

    /**
     * @brief       compute winnowed minimizers from a given sequence and add to the index
     * @param[out]  minimizerIndex  minimizer table storing minimizers and their position as we compute them
     * @param[in]   seq             pointer to input sequence
     * @param[in]   len             length of input sequence
     * @param[in]   kmerSize
     * @param[in]   windowSize
     * @param[in]   seqCounter      current sequence number, used while saving the position of minimizer
     */
    template <typename T>
      inline void addMinimizers(std::vector<T> &minimizerIndex, 
          char* seq, offset_t len,
          int kmerSize, 
          int windowSize,
          int alphabetSize,
          seqno_t seqCounter)
      {
        /**
         * Double-ended queue (saves minimum at front end)
         * Saves pair of the minimizer and the position of hashed kmer in the sequence
         * Position of kmer is required to discard kmers that fall out of current window
         */
        std::deque< std::pair<MinimizerInfo, offset_t> > Q;

        makeUpperCase(seq, len);

        //Compute reverse complement of seq
        char* seqRev = new char[len];

        if(alphabetSize == 4) //not protein
          CommonFunc::reverseComplement(seq, seqRev, len);

        for(offset_t i = 0; i < len - kmerSize + 1; i++)
        {
          //The serial number of current sliding window
          //First valid window appears when i = windowSize - 1
          offset_t currentWindowId = i - windowSize + 1;

          //Hash kmers
          hash_t hashFwd = CommonFunc::getHash(seq + i, kmerSize); 
          hash_t hashBwd;

          if(alphabetSize == 4)
            hashBwd = CommonFunc::getHash(seqRev + len - i - kmerSize, kmerSize);
          else  //proteins
            hashBwd = std::numeric_limits<hash_t>::max();   //Pick a dummy high value so that it is ignored later

          //Consider non-symmetric kmers only
          if(hashBwd != hashFwd)
          {
            //Take minimum value of kmer and its reverse complement
            hash_t currentKmer = std::min(hashFwd, hashBwd);

            //Check the strand of this minimizer hash value
            auto currentStrand = hashFwd < hashBwd ? strnd::FWD : strnd::REV;

            //If front minimum is not in the current window, remove it
            while(!Q.empty() && Q.front().second <=  i - windowSize)
              Q.pop_front();

            //Hashes less than equal to currentKmer are not required
            //Remove them from Q (back)
            while(!Q.empty() && Q.back().first.hash >= currentKmer) 
              Q.pop_back();

            //Push currentKmer and position to back of the queue
            //-1 indicates the dummy window # (will be updated later)
            Q.push_back( std::make_pair(
                  MinimizerInfo{currentKmer, seqCounter, -1, currentStrand},
                  i)); 

            //Select the minimizer from Q and put into index
            if(currentWindowId >= 0)
            {
              //We save the minimizer if we are seeing it for first time
              if(minimizerIndex.empty() || minimizerIndex.back() != Q.front().first)
              {
                //Update the window position in this minimizer
                //This step also ensures we don't re-insert the same minimizer again
                Q.front().first.wpos = currentWindowId;     
                minimizerIndex.push_back(Q.front().first);
              }
            }
          }
        }

#ifdef DEBUG
        //std::cout << "INFO, skch::CommonFunc::addMinimizers, inserted minimizers for sequence id = " << seqCounter << "\n";
#endif

        delete [] seqRev;
      }

   /**
     * @brief           Functor for comparing tuples by single index layer
     * @tparam layer    Tuple's index which is used for comparison
     * @tparam op       comparator, default as std::less
     */
    template <size_t layer, template<typename> class op = std::less>
      struct TpleComp
      {
        //Compare two tuples using their values
        template<typename T>
          bool operator() (T const &t1, T const &t2)
          {
            return op<typename std::tuple_element<layer, T>::type>() (std::get<layer>(t1), std::get<layer>(t2));
          }
      };

    /**
     * @brief                   computes the total size of reference in bytes
     * @param[in] refSequences  vector of reference files
     * @return                  total size
     */
    inline uint64_t getReferenceSize(const std::vector<std::string> &refSequences)
    {
      uint64_t count = 0;

      for(auto &f : refSequences)
      {
        //Open the file as binary, and set the position to end
        std::ifstream in(f, std::ifstream::ate | std::ifstream::binary);

        //the position of the current character
        count += (uint64_t)(in.tellg());
      }

      return count;
    }

  }
}

#endif
