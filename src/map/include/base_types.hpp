/**
 * @file    base_types.hpp
 * @brief   Critical type defintions for mapping algorithm
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef BASE_TYPES_MAP_HPP 
#define BASE_TYPES_MAP_HPP

#include <tuple>
#include <vector>

namespace skch
{
  typedef uint32_t hash_t;    //hash type
  typedef int offset_t;       //position within sequence
  typedef int seqno_t;        //sequence counter in file
  typedef int8_t strand_t;   //sequence strand 
  typedef int8_t side_t;   //sequence strand 

  //C++ timer
  typedef std::chrono::high_resolution_clock Time;

  //Information about each minimizer
  struct MashimizerInfo
  {
    hash_t hash;                              //hash value
    seqno_t seqId;                            //sequence or contig id
    offset_t wpos;                            //First (left-most) window position when the minimizer is saved
    offset_t wpos_end;
    strand_t strand;                          //strand information

    //Lexographical equality comparison
    bool operator ==(const MashimizerInfo& x) const {
      return std::tie(hash, seqId, wpos, strand) 
        == std::tie(x.hash, x.seqId, x.wpos, x.strand);
    }

    bool operator !=(const MashimizerInfo& x) const {
      return std::tie(hash, seqId, wpos, strand) 
        != std::tie(x.hash, x.seqId, x.wpos, x.strand);
    }

    static bool equalityByHash(const MashimizerInfo& x, const MashimizerInfo& y) {
      return x.hash == y.hash;
    }

    static bool lessByHash(const MashimizerInfo& x, const MashimizerInfo& y) {
      return x.hash < y.hash;
    }

    // Sort based on start point
    bool operator <(const MashimizerInfo& x) const {
      return std::tie(seqId, wpos) 
        < std::tie(x.seqId, x.wpos);
    }
  };


  //Type for map value type used for
  //L1 stage lookup index
  struct MinimizerMetaData
  {
    seqno_t seqId;          //sequence or contig id
    offset_t wpos;          //window position (left-most window)
    strand_t strand;        //strand information

    bool operator <(const MinimizerMetaData& x) const {
      return std::tie(seqId, wpos, strand) 
        < std::tie(x.seqId, x.wpos, x.strand);
    }
  };

  template <class ostream>
  ostream& operator<< (ostream& os, const MashimizerInfo& mi) {
    os << mi.seqId <<  "\t" << std::to_string(mi.strand) << "\t" << mi.wpos << "\t" << mi.wpos_end << "\t" << mi.hash;
    return os;
  };

  // Enum for tracking which side of an interval a point represents
  enum side : side_t
  {
    OPEN = 1,  
    CLOSE = -1
  };  

  // Endpoints for mashimizer intervals
  struct IntervalPoint
  {
    side_t side;
    seqno_t seqId;
    offset_t pos;
    strand_t strand;
    hash_t hash;

    // Sort interval points. 
    // For a pair of points at the same seqId/pos, the end point should be first
    bool operator <(const IntervalPoint& x) const {
      return std::tie(seqId, pos, side) 
        < std::tie(x.seqId, x.pos, x.side);
    }
  };

  template <class ostream>
  ostream& operator<< (ostream& os, const IntervalPoint& mi) {
    os << mi.seqId <<  "\t" << std::to_string(mi.strand) << "\t" << mi.pos << "\t" << (mi.side == side::OPEN ? "OPEN" : "CLOSE");
    return os;
  };


  typedef hash_t MinimizerMapKeyType;
  typedef std::vector<MashimizerInfo> MinimizerMapValueType;

  //Metadata recording for contigs in the reference DB
  struct ContigInfo
  {
    std::string name;       //Name of the sequence
    offset_t len;           //Length of the sequence
  };

  //Label tags for strand information
  enum strnd : strand_t
  {
    FWD = 1,  
    AMBIG = 0,
    REV = -1
  };  

  enum event : int
  {
    BEGIN = 1,
    END = 2
  };

  //filter mode in mashmap
  enum filter : int
  {
    MAP = 1,                              //filter by query axis
    ONETOONE = 2,                         //filter by query axis and reference axis
    NONE = 3                              //no filtering
  };

  //Fragment mapping result
  //Do not save variable sized objects in this struct
  struct MappingResult
  {
    offset_t queryLen;                                  //length of the query sequence
    offset_t refStartPos;                               //start position of the mapping on reference
    offset_t refEndPos;                                 //end pos
    offset_t queryStartPos;                             //start position of the query for this mapping
    offset_t queryEndPos;                               //end position of the query for this mapping
    seqno_t refSeqId;                                   //internal sequence id of the reference contig
    seqno_t querySeqId;                                 //internal sequence id of the query sequence
    float nucIdentity;                                  //calculated identity
    float nucIdentityUpperBound;                        //upper bound on identity (90% C.I.)
    int sketchSize;                                     //sketch size
    int conservedSketches;                              //count of conserved sketches
    strand_t strand;                                    //strand

                                                        //--for split read mapping

    offset_t splitMappingId;                            // To identify split mappings that are chained
    int discard;                                        // set to 1 for deletion

    offset_t qlen() {                                   //length of this mapping on query axis 
      return queryEndPos - queryStartPos + 1;
    }

    offset_t rlen() {                                   //length of this mapping on reference axis 
      return refEndPos - refStartPos + 1;
    }
  };

  typedef std::vector<MappingResult> MappingResultsVector_t;

  //Container to save copy of kseq object
  struct InputSeqContainer
  {
    seqno_t seqCounter;                 //sequence counter
    offset_t len;                       //sequence length             
    std::string seq;                    //sequence string
    std::string seqName;                //sequence id

    /*
     * @brief               constructor
     * @param[in] kseq_seq  complete read or reference sequence
     * @param[in] kseq_id   sequence id name
     * @param[in] len       length of sequence
     */
    InputSeqContainer(const char * kseq_seq_, const char * kseq_id_, offset_t len_, seqno_t seqcount_)
    {
      this->seq = std::string{kseq_seq_, std::size_t(len_)};
      this->seqName = std::string{kseq_id_};
      this->len = len_;
      this->seqCounter = seqcount_; 
    }
  };

  //Output type of map function
  struct MapModuleOutput
  {
    MappingResultsVector_t readMappings;  //read mapping coordinates
    std::string qseqName;                 //query sequence id

    //Function to erase all output mappings
    void reset()
    {
      this->readMappings.clear();
    }
  };

  //Information about fragment sequence during L1/L2 mapping
  template <typename MinimizerVec>
    struct QueryMetaData
    {
      char *seq;                          //query sequence pointer 
      seqno_t seqCounter;                 //query sequence counter
      std::string seqName;
      offset_t len;                       //length of this query sequence
      offset_t startPos;
      int sketchSize;                     //sketch size
      MinimizerVec minimizerTableQuery;   //Vector of minimizers in the query 
      MinimizerVec seedHits;              //Vector of seed hits in the reference
    };
}

#endif
