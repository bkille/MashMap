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

//External includes
#include "common/murmur3.h"
#include "common/prettyprint.hpp"
#include "map/include/base_types.hpp"

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
        assert(sketched_vals.size() <= sketchSize);
        mashimizerIndex.reserve(sketched_vals.size());
        std::transform(sketched_vals.begin(), sketched_vals.end(), std::back_inserter(mashimizerIndex), 
            [](auto& pair) {pair.second.strand = pair.second.strand >= 0 ? strnd::FWD : strnd::REV; return pair.second;});
        //std::for_each(mashimizerIndex.begin(), mashimizerIndex.end(), [](auto& mi) {mi.strand = (mi.strand >= 0 ? strnd::FWD : strnd::REV);});
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


        std::deque< std::pair<hash_t, offset_t> > Q;
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
          offset_t currentWindowId = i - windowSize + 1;

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

          //Consider non-symmetric kmers only
          if(hashBwd != hashFwd)
          {
            //Take minimum value of kmer and its reverse complement
            hash_t currentKmer = std::min(hashFwd, hashBwd);

            //Check the strand of this minimizer hash value
            auto currentStrand = hashFwd < hashBwd ? strnd::FWD : strnd::REV;

            //If front minimum is not in the current window, remove it
            if (!Q.empty() && Q.front().second <=  i - windowSize) {
              const hash_t leaving_hash = Q.front().first;

              // If the hash that is getting popped off is still in the window and it is now leaving the window 
              // wpos != -1 and wpos_end == -1 --> still in window
              if (sortedWindow[leaving_hash].first.wpos != -1 and sortedWindow[leaving_hash].first.wpos_end == -1 && sortedWindow[leaving_hash].second == 1) {
                sortedWindow[leaving_hash].first.wpos_end = currentWindowId;
                mashimizerIndex.push_back(sortedWindow[leaving_hash].first);
              }

              // Remove hash
              // Popping back for now, to retain the window. Should probably change
              // the value type to be a pair of MI and count
              sortedWindow[leaving_hash].second -= 1;
              if (sortedWindow[leaving_hash].second == 0) {
                if (leaving_hash == piv.p->first) {
                  std::advance(piv.p, -1);
                  piv.rank -= 1;
                }
                else if (leaving_hash < piv.p->first || sortedWindow.size() < sketchSize) {
                  piv.rank -= 1;
                }
                sortedWindow.erase(leaving_hash);
              }
              Q.pop_front();
            }

            // Add current hash to window
            if (sortedWindow.size() < sketchSize*2+20 || currentKmer <= std::prev(sortedWindow.end())->first) {
                Q.push_back(std::make_pair(currentKmer, i)); 
                if (sortedWindow[currentKmer].second == 0) {
                    auto mi = MashimizerInfo{currentKmer, seqCounter, -1, -1, currentStrand};
                    sortedWindow[currentKmer].first = mi;
                }
                sortedWindow[currentKmer].second += 1;
            }


            //Select the minimizer from Q and put into index
            if(currentWindowId >= 0)
            {
              bool insert_mi = mashimizerIndex.empty();
              // There are two cases in which we add a new mashimizer to the index
              // (1) currentKmer <= piv.p->first || sortedWindow.size() <= sketchSize && currentKmer is new
              // (2) We just removed an element from the sketch, and therefore shifted a new element
              // in. 
              if (sortedWindow.size() <= sketchSize || currentKmer < piv.p->first)
              {
                // New kmer in sketch
                if (sortedWindow[currentKmer].second == 1) {
                  if (sortedWindow.size() <= sketchSize) {
                      piv.p = std::prev(sortedWindow.end());
                      piv.rank = sortedWindow.size();
                  } else { 
                      if (piv.rank == sketchSize)
                        std::advance(piv.p, -1);
                      piv.rank = sketchSize;
                  }
                }
                //Update the window position in this mashimizer
                sortedWindow[currentKmer].first.wpos = sortedWindow[currentKmer].first.wpos >= 0 ?
                    sortedWindow[currentKmer].first.wpos : currentWindowId;     
              } 
              if (sortedWindow.size() >= sketchSize && piv.rank == sketchSize - 1) {
                std::advance(piv.p, 1);
                piv.rank += 1;
                // hash on the border just got added back in. In some cases, it could've been
                // in from the start of this loop iteration but gotten kicked out by the new kmer
                piv.p->second.first.wpos = piv.p->second.first.wpos >= 0 ? piv.p->second.first.wpos : currentWindowId;     
              } else {}
                // Check if we just removed a hash from the sketch
              if (sortedWindow.size() > sketchSize) {
                auto& border_mi_it = std::next(piv.p)->second.first;
                if (border_mi_it.wpos != -1 && border_mi_it.wpos_end == -1) {
                  border_mi_it.wpos_end = currentWindowId;
                  mashimizerIndex.push_back(
                          MashimizerInfo{border_mi_it.hash, border_mi_it.seqId, border_mi_it.wpos, border_mi_it.wpos_end, border_mi_it.strand});
                  // resest mi info
                  border_mi_it.wpos = -1;
                  border_mi_it.wpos_end = -1;
                }
              }
            } else {
                if (sortedWindow[currentKmer].second == 1) {
                    // Seeing kmer for the first time
                    if (sortedWindow.size() <= sketchSize) {
                        piv.p = std::prev(sortedWindow.end());
                        piv.rank = sortedWindow.size();
                    } else if (currentKmer < piv.p->first) {
                        std::advance(piv.p, -1);
                    }
                }
            }

          }
//#ifdef DEBUG
          //if (i % 10000000 == 0 and i != 0) {
              //std::cout << i << std::endl;
              //std::cout << piv.rank << ", " << sortedWindow.size() << ", " << i << std::endl;
          //}
          //if ((sortedWindow.size() > 0 ? std::distance(sortedWindow.begin(), piv.p) + 1 : 0) != piv.rank) {
              //std::cout << "Actual rank = " 
                  //<< (sortedWindow.size() > 0 ? std::distance(sortedWindow.begin(), piv.p) + 1 : 0 )
                  //<< "\tAnnotated rank = " << piv.rank << std::endl;
              //exit(1);
          //}
          //if (piv.rank != std::min<uint64_t>(sortedWindow.size(), sketchSize) ){
              //std::cout << "Actual rank = " 
                  //<< (sortedWindow.size() > 0 ? std::distance(sortedWindow.begin(), piv.p) + 1 : 0 )
                  //<< "\tAnnotated rank = " << piv.rank << std::endl;
              //exit(1);
          //}
//#endif
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

#ifdef DEBUG
        std::sort(mashimizerIndex.begin(), mashimizerIndex.end(), [](auto& l, auto& r) {
            if (l.hash != r.hash) return l.hash < r.hash; else return l.wpos < r.wpos;});
        auto prev_hash = 0;
        for (auto iter = mashimizerIndex.begin() + 1; iter != mashimizerIndex.end(); iter++) {
          if (iter->hash != (iter-1)->hash) {
            continue;
          }
          assert((iter-1)->wpos_end <= iter->wpos);
        }
#endif 

        // TODO is this necessary?
        std::sort(mashimizerIndex.begin(), mashimizerIndex.end(), [](auto& l, auto& r) {return l.wpos < r.wpos;});
        mashimizerIndex.erase(std::unique(mashimizerIndex.begin(), mashimizerIndex.end(), 
                    [](auto& l, auto& r) {
                        return (l.wpos == r.wpos) && (l.hash == r.hash);
            }),
            mashimizerIndex.end());
#ifdef DEBUG
        std::cout << "INFO, skch::CommonFunc::addMinimizers, inserted minimizers for sequence id = " << seqCounter << "\n";
        std::cout << "INFO, skch::CommonFunc::addMinimizers, length of sequence  = " << len << "\n";
        assert(std::all_of(mashimizerIndex.begin(), mashimizerIndex.end(), [](auto& mi) {return mi.wpos >= 0;}));
        assert(std::all_of(mashimizerIndex.begin(), mashimizerIndex.end(), [](auto& mi) {return mi.wpos_end >= 0;}));
        std::vector<MashimizerInfo> endpos_heap;
        auto heap_cmp = [](auto& l, auto& r) {return l.wpos_end >= r.wpos_end;};
        for (auto& mi : mashimizerIndex) {
          while (!endpos_heap.empty() && endpos_heap.front().wpos_end <= mi.wpos) {
            std::pop_heap(endpos_heap.begin(), endpos_heap.end(), heap_cmp); 
            endpos_heap.pop_back();
          }
          endpos_heap.push_back(mi);
          std::push_heap(endpos_heap.begin(), endpos_heap.end(), heap_cmp);
          assert(endpos_heap.size() <= sketchSize);
        }
#endif
        //std::cout << "BEFORE " << mashimizerIndex.size() << "\n";
        //std::cout << "AFTER " << mashimizerIndex.size() << "\n";
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
        std::cout << "INFO, skch::CommonFunc::addMinimizers, inserted minimizers for sequence id = " << seqCounter << "\n";
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
