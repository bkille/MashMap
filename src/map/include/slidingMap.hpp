/**
 * @file    slidingMap.hpp
 * @brief   implements ordered map to compute Jaccard
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef SLIDING_MAP_HPP 
#define SLIDING_MAP_HPP

#include <vector>
#include <algorithm>
#include <map>

//Own includes
#include "map/include/base_types.hpp"

//External includes
//#include "assert.hpp"

namespace skch
{
  /**
   * @class     skch::SlideMapper
   * @brief     L1 and L2 mapping stages
   */
  template <typename Q_Info>
    class SlideMapper
    {

      private:

        //Metadata for the minimizers saved in sliding ordered map during L2 stage
        struct slidingMapContainerValueType
        {
          unsigned int freq;
          strand_t q_strand;
          strand_t r_strand;
          bool is_query_mashimizer;
        };

        //Container type for saving read sketches during L1 and L2 both
        typedef Sketch::MI_Type MinVec_Type;

        typedef Sketch::MIIter_t MIIter_t;

        //reference to query's metadata
        const Q_Info &Q;

        //Define a Not available position marker
        static const offset_t NAPos = std::numeric_limits<offset_t>::max();

        //Ordered map to save unique sketch elements, and associated value as 
        //a pair of its occurrence in the query and the reference
        typedef std::map< hash_t, slidingMapContainerValueType> MapType;

        //Iterator pointing to the smallest 's'th element in the map
        typename MapType::iterator pivot;

        //Label status while inserting reference minimizer
        enum IN : int
        {
          //reference minimizer is inserted into map as a new map entry
          UNIQ = 1,

          //reference minimizer is coupled with a query minimizer, 
          //previously, there was no ref. minimizer at this entry
          CPLD = 2, 

          //ref. minimizer just revises the hash position of already 
          //existing reference minimizer
          REV = 3
        };  

        //Label status while deleting reference minimizer
        enum OUT : int
        {
          //entry in the map is deleted
          DEL = 1,

          //just the reference minimizer is updated to null
          UPD = 2,

          //Nothing changed in the map
          NOOP = 3
        };

      public:

        MapType slidingWindowMinhashes;

        //Count of shared sketch elements between query and the reference
        //Updated after insert or delete operation on map 
        int sharedSketchElements;
        int strand_votes;


        //Delete default constructor
        SlideMapper() = delete;

        /**
         * @brief                 constructor
         * @param[in]   Q         query meta data
         */
        SlideMapper(Q_Info &Q_) :
          Q(Q_),
          pivot(this->slidingWindowMinhashes.end()),
          sharedSketchElements(0),
          strand_votes(0)
        {
          this->init();
        }

      private:

        /**
         * @brief       Fills map with minimum 's' minimizers in the query
         */
        inline void init()
        {
          //Insert query sketch elements to map
          for(auto it = Q.minimizerTableQuery.begin(); it != Q.minimizerTableQuery.end(); it++)
          {
            // TODO I don't think the Q.minimizerTableQuery vec is sorted...
            this->slidingWindowMinhashes.emplace_hint(slidingWindowMinhashes.end(), it->hash, slidingMapContainerValueType {0, it->strand, 0, true});

            //slidingWindowMinhashes[it->hash] = slidingMapContainerValueType {1, true};
          }

          //std::cout << slidingWindowMinhashes.size() << " min hashes to start\n";
          //Point pivot to last element in the map
          this->pivot = std::prev(this->slidingWindowMinhashes.end());

          //Current count of shared sketch elements is zero
          this->sharedSketchElements = 0;
        }

      public:

        /**
         * @brief               insert a minimizer from the reference sequence into the map
         * @param[in]   m       reference minimizer to insert
         */
        inline void insert_mashimizer(const skch::MashimizerInfo& mi)
        {
          const hash_t hashVal = mi.hash;
          int status;

          //if hash doesn't exist in the map, add to it
          if(slidingWindowMinhashes.find(hashVal) == slidingWindowMinhashes.end())
          {
            slidingWindowMinhashes[hashVal] = slidingMapContainerValueType {1, 0, 0, false};   //add the hash to window
            status = IN::UNIQ;
            //std::cout << "Unique insert " << mi.hash << " @ " << mi.wpos << std::endl; 
          }
          else
          {
            status = (slidingWindowMinhashes[hashVal].is_query_mashimizer && slidingWindowMinhashes[hashVal].freq == 0) ? IN::CPLD 
              : IN::REV;
            if (status == IN::CPLD) {
              //std::cout << "Pairing " << mi.hash << " @ " << mi.wpos << std::endl; 
            }

            //if hash already exists in the map, just revise it
            slidingWindowMinhashes[hashVal].freq += 1;
          }

          updateCountersAfterInsert(status, mi);

          assert(std::distance(slidingWindowMinhashes.begin(), pivot) == Q.sketchSize - 1);
          assert(this->sharedSketchElements >= 0);
          assert(this->sharedSketchElements <= Q.sketchSize);
        }

        /**
         * @brief               delete a minimizer from the reference sequence from the map
         * @param[in]   m       reference minimizer to remove
         */
        void delete_mashimizer(const skch::MashimizerInfo& mi)
        {
          int status;
          const hash_t hashVal = mi.hash;
          bool pivotDeleteCase = false;
          
          //assert(this->slidingWindowMinhashes.find(hashVal) != this->slidingWindowMinhashes.end(), "Can't find hash to delete", hashVal);
          // End point may not have had an open point
          if (this->slidingWindowMinhashes.find(hashVal) == this->slidingWindowMinhashes.end()) {
              return;
          }

          this->slidingWindowMinhashes[hashVal].freq -= 1;
          if (this->slidingWindowMinhashes[hashVal].freq > 0) {
              return;
          }

          if(!this->slidingWindowMinhashes[hashVal].is_query_mashimizer)
          {
            //Handle pivot deletion as a separate case
            if(hashVal == pivot->first)
            {
              pivot++;
              if(pivot->second.is_query_mashimizer && pivot->second.freq > 0)
              {
                //std::cout << "Special pivot case\n";
                if (this->pivot->second.r_strand == 0) {
                  //std::cout << "ERROR\n";
                }

                this->sharedSketchElements += 1;
                this->strand_votes += this->pivot->second.r_strand * this->pivot->second.q_strand;
              }
              pivotDeleteCase = true;
            }

            //std::cout << "Removing " << mi.hash <<  " @ " << mi.wpos << std::endl; 
            this->slidingWindowMinhashes.erase(hashVal);              //Remove the entry from the map
            status = OUT::DEL; 
          }
          else
          {
            //std::cout << "Unpairing " << mi.hash << " @ " << mi.wpos << std::endl; 
            status = OUT::UPD; 
          }

          if(!pivotDeleteCase) 
            updateCountersAfterDelete(status, mi);

          assert(std::distance(slidingWindowMinhashes.begin(), pivot) == Q.sketchSize - 1);
          assert(this->sharedSketchElements >= 0);
          assert(this->sharedSketchElements <= Q.sketchSize);
        }


      private:

        /**
         * @brief             logic to update internal counters after insert to map
         * @param[in] status  insert status
         * @param[in] m       reference minimizer that was inserted
         */
        void updateCountersAfterInsert(int status, const skch::MashimizerInfo& mi)
        {
          hash_t hash = mi.hash;
          //Revise internal counters
          if(hash <= this->pivot->first)
          {
            if(status == IN::CPLD)
            {
              //Increase count of shared sketch elements by 1
              this->sharedSketchElements += 1;
              this->slidingWindowMinhashes[hash].r_strand = mi.strand;
              if (this->slidingWindowMinhashes[hash].r_strand == 0) {
                //std::cout << hash <<  " ERROR3\n";
              }
              this->strand_votes += this->slidingWindowMinhashes[hash].q_strand * mi.strand;
            }
            else if(status == IN::UNIQ)
            {
              if (mi.hash < pivot->first) {
                if(this->pivot->second.is_query_mashimizer && this->pivot->second.freq > 0) {
                  this->sharedSketchElements -= 1;
                  this->strand_votes -= this->pivot->second.r_strand * this->pivot->second.q_strand;
                }


                //Pivot needs to be decremented
                std::advance(this->pivot, -1);
              }
            }
            else if(status == IN::REV)
            {
              //Do nothing
            }
          } else {
            if(status == IN::CPLD)
            {
              this->slidingWindowMinhashes[hash].r_strand = mi.strand;
            }
          }
        }

        /**
         * @brief             logic to update internal counters after delete to map
         * @param[in] status  insert status
         * @param[in] m       reference minimizer that was inserted
         */
        void updateCountersAfterDelete(int status, const skch::MashimizerInfo& mi)
        {
          //Revise internal counters
          if(mi.hash <= this->pivot->first)
          {
            if(status == OUT::UPD)
            {
              //Decrease count of shared sketch elements by 1
              this->sharedSketchElements -= 1;
              this->strand_votes -= this->slidingWindowMinhashes[mi.hash].q_strand * this->slidingWindowMinhashes[mi.hash].r_strand;
              this->slidingWindowMinhashes[mi.hash].r_strand = 0;
            }
            else if(status == OUT::DEL)
            {
              //Pivot needs to be advanced
              std::advance(this->pivot, 1);

              if(this->pivot->second.is_query_mashimizer && this->pivot->second.freq > 0)
              {
                this->sharedSketchElements += 1;
                this->strand_votes += pivot->second.q_strand * pivot->second.r_strand;
                if (this->pivot->second.r_strand == 0) {
                  //std::cout << pivot->first <<  "ERROR2\n";
                }
              }
            }
            else if(status == OUT::NOOP)
            {
              //Do nothing
            }
          }
        }

    };
}

#endif
