#include <vector>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/AAIndex.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>

using namespace OpenMS;
using namespace std;

class inner_Node
{

  private:

  Peak1D fragment_;

  const AASequence* peptide_ptr_;

  double prec_weight_;

  double frag_weight_;

  public:

  inner_Node() : peptide_ptr_(NULL), fragment_({}), prec_weight_(0), frag_weight_(0){}

  inner_Node(const Peak1D& fragment, const AASequence* peptide) : peptide_ptr_(peptide), fragment_(fragment), prec_weight_((*peptide).getMonoWeight()), frag_weight_(fragment.getMZ()){}

  inner_Node(const inner_Node& n)
  {

    fragment_ = n.fragment_;

    peptide_ptr_ = n.peptide_ptr_;

    prec_weight_ = n.prec_weight_;

    frag_weight_ = n.frag_weight_;

  }

  inner_Node(double prec_weight) : fragment_({}), peptide_ptr_(NULL), prec_weight_(prec_weight), frag_weight_(0){}

  bool operator< (const inner_Node& r) const
  {
    
    if (prec_weight_ == r.prec_weight_)
    {

      return(frag_weight_ < r.frag_weight_);

    }
    else
    {

      return (prec_weight_ < r.prec_weight_);

    }    

  }

  double getfragWeight()const
  {

    return frag_weight_;

  }

  double getprecWeight()const
  {

    return prec_weight_;

  }

  const AASequence& getPrecursor() const
  {

    return *peptide_ptr_;

  }

};

class outer_Node
{

  private:

  std::set<inner_Node> nodes_;

  size_t size_;

  double min_;

  public:

  outer_Node()
  {

    nodes_ = {};

    size_ = 0;

    min_ = 0;

  }

  outer_Node(double min) : nodes_({}), size_(0), min_(min){}

  outer_Node(const outer_Node& n)
  {

    nodes_ = n.nodes_;

    size_ = n.size_;

    min_ = n.min_;

  }

  size_t getSize() const
  {

    return nodes_.size();

  }

  void insert(const inner_Node& n)
  {

    nodes_.insert(n);

    if(size_ == 0){

      min_ = n.getfragWeight();

    }

    size_++;

  }

  void clear()
  {

    nodes_ = {};

    size_ = 0;

    min_ = 0;

  }

  bool operator<(const outer_Node& r) const
  {

    return (min_ < r.min_);

  }

  double getMin() const
  {

    return min_;

  }

  const std::set<inner_Node>& getNodes() const
  {

    return nodes_;

  }

};

class Tree
{

  private:

  std::set<outer_Node> nodes_;

  std::vector<AASequence> all_peptides_;

  size_t size_;

  struct Temp_frag
  {
    const AASequence* prec_ptr_;


    Peak1D frag_;


    Temp_frag(const AASequence& prec, const Peak1D& frag) : prec_ptr_(&prec), frag_(frag){}


    bool operator< (const Temp_frag& r) const
    {
      return (frag_.getMZ() < r.frag_.getMZ());
    }
  };


  std::vector<AASequence> load_and_digest(std::string fastafile, std::string digestor_enzyme)
  {
    //loading
    FASTAFile file;
    std::vector<FASTAFile::FASTAEntry> entries;
    file.load(fastafile, entries);


    //setting digestor
    ProteaseDigestion digestor;
    digestor.setEnzyme(digestor_enzyme);
    std::vector<AASequence> peptides;

    std::vector<AASequence> all_peptides{};


    for (auto entry : entries)
    {
    //digesting
    digestor.digest(AASequence::fromString(entry.sequence), peptides);


    all_peptides.insert(all_peptides.end(), std::make_move_iterator(peptides.begin()), std::make_move_iterator(peptides.end())); // Move: make_move_iterator
    }
    return all_peptides;
  }

  void frag_merge(const vector<int>& chunks, vector<Temp_frag>& input)
  {
    for (size_t i=0; i<chunks.size();i++)
    {

    }
    
  }

  const std::vector<Temp_frag> generate_frag()
  {
    TheoreticalSpectrumGenerator tsg;
    PeakSpectrum b_y_ions;
    std::vector<Temp_frag> all_frags;
    std::vector<int> chunk_start = {0};

    // Bottleneck
    for(const auto& pep : all_peptides_)
    {
      
      if(pep.toString().find('X') == std::string::npos)
      {
        
        tsg.getSpectrum(b_y_ions, pep, 1, 1);      

      for(const auto& frag : b_y_ions)
      {
        all_frags.emplace_back(pep, frag); 
      }
      chunk_start.emplace_back(all_frags.size()-1);
      b_y_ions.clear(true);
      }
    }
    
    frag_merge(chunk_start, all_frags);

    return all_frags;
  }

  public:

  Tree()
  {
    nodes_ = {};
    size_ = 0;
  }

  Tree(std::string fastafile, std::string digestor_enzyme, size_t bucketsize){


    all_peptides_ = load_and_digest(fastafile, digestor_enzyme);
    
    const std::vector<Temp_frag> all_frags = generate_frag();
    
    size_t i = 0;
    size_ = 0;
    outer_Node new_node{};

    for(const auto& frag : all_frags)
    {

      if(i < (bucketsize-1))
      {
        
        inner_Node new_inner_node{frag.frag_, frag.prec_ptr_};

        new_node.outer_Node::insert(new_inner_node);
        i++;

      }
      else
      {

        inner_Node new_inner_node{frag.frag_, frag.prec_ptr_};

        new_node.outer_Node::insert(new_inner_node);
        
        nodes_.insert(new_node);
        new_node.outer_Node::clear();
        i = 0;
        size_++;

      }

    }
    nodes_.insert(new_node);
    size_++;

  }
  
  size_t getSize() const
  {

    return size_;

  }

  void search(const MSSpectrum& spectrum, const double tolerance, vector<AASequence>& candidates) const 
  {

    candidates.clear();

    std::vector<Precursor> Precursor = spectrum.getPrecursors();
    double prec_mz = 0;
    if(Precursor.size() == 1)
    {

      prec_mz = Precursor[0].getUnchargedMass();      

    }
    else{}

    inner_Node search_inner_upper{(prec_mz + tolerance)};
    inner_Node search_inner_lower{(prec_mz - tolerance)};

    for(auto peak : spectrum)
    {
      
      outer_Node search_outer_upper{(peak.getMZ() + tolerance)};
      outer_Node search_outer_lower{(peak.getMZ() - tolerance)};

      auto outer_upper = nodes_.upper_bound(search_outer_upper);
      auto outer_lower = nodes_.lower_bound(search_outer_lower);

      for(auto it = outer_lower; it != outer_upper; it++)
      {

        auto inner_upper = (*it).getNodes().upper_bound(search_inner_upper);
        auto inner_lower = (*it).getNodes().lower_bound(search_inner_lower);

        for(auto it2 = inner_lower; it2 != inner_upper; it2++)
        {

          candidates.push_back((*it2).getPrecursor());

        }

      }

    }

  }

};

int main(){

  Tree test("/buffer/ag_bsc/pmsb_23/heike/fasta_testing/philosopher/test/db/uniprot/hsa-reviewed-2019-02-04.fasta", "Trypsin", 16);
  // Tree test("/buffer/ag_bsc/pmsb_23/max_alcer/sage/tests/Q99536.fasta", "Trypsin", 16);
  vector<AASequence> candidates;

  MzMLFile f;

  MSExperiment exp;

  f.load("/buffer/ag_bsc/pmsb_23/max_alcer/sage/tests/LQSRPAAPPAPGPGQLTLR.mzML", exp);

  test.search(*exp.begin(), 0.01, candidates);

  //Tree test("/buffer/ag_bsc/pmsb_23/max_alcer/OpenMS/src/tests/class_tests/openms/data/Sequest_test2.fasta", "Trypsin", 16);
/*
  outer_Node test_o = test.test_return();
  cout << "\n";
  for(auto i : test_o.getNodes()){
    
    cout << i.getprecWeight() << "\n";

  } */

  // std::cout << candidates[0].toString() << " " << candidates[1].toString() << std::endl;
  return 0;

}

