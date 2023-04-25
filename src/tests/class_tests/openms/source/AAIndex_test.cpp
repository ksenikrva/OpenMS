#include <vector>
#include <omp.h>
#include <OpenMS/ANALYSIS/ID/SearchDatabase.h>
#include <OpenMS/ANALYSIS/ID/SimpleSearchEngineAlgorithm.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
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
/*
class inner_Node
{

  private:

  const AASequence* peptide_ptr_{};

  Peak1D fragment_{};  

  double prec_weight_{};

  double frag_weight_{};

  public:

  inner_Node() = delete;

  inner_Node(const AASequence* peptide, const Peak1D& fragment) : peptide_ptr_(peptide), fragment_(fragment), prec_weight_((*peptide).getMonoWeight()), frag_weight_(fragment.getMZ()){}

  inner_Node(const inner_Node& n) = default; 

  inner_Node(double prec_weight) : prec_weight_(prec_weight){}

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

  std::set<inner_Node> nodes_{};

  size_t size_{};

  double min_{};

  public:

  outer_Node() = default;

  outer_Node(double min) : min_(min){}

  outer_Node(const outer_Node& n) = default;

  size_t getSize() const
  {

    return size_;

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

    *this = {};

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

  const std::set<Temp_frag> generate_frag()
  {
    TheoreticalSpectrumGenerator tsg;
    PeakSpectrum b_y_ions;
    std::set<Temp_frag> all_frags;


    // Bottleneck
    for(const auto& pep : all_peptides_)
    {
      
      if(pep.toString().find('X') == std::string::npos)
      {
        
        tsg.getSpectrum(b_y_ions, pep, 1, 1);      

      for(const auto& frag : b_y_ions)
      {
        all_frags.insert({pep, frag}); 
      }
      b_y_ions.clear(true);
      }
    }
    return all_frags;
  }

  public:

  Tree() = delete;

  Tree(std::string fastafile, std::string digestor_enzyme, size_t bucketsize){


    all_peptides_ = load_and_digest(fastafile, digestor_enzyme);
    
    const std::set<Temp_frag> all_frags = generate_frag();
    
    size_t i = 0;
    size_ = 0;
    outer_Node new_node{};

    for(const auto& frag : all_frags)
    {

      if(i < (bucketsize-1))
      {
        
        inner_Node new_inner_node{frag.prec_ptr_, frag.frag_};

        new_node.outer_Node::insert(new_inner_node);
        i++;

      }
      else
      {

        inner_Node new_inner_node{frag.prec_ptr_, frag.frag_};

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
*/

// OpenMS::SearchDatabase sdb_test(std::vector<OpenMS::FASTAFile::FASTAEntry>{});

void fragment_merge(int first, int last, const std::vector<int>& chunks, std::vector<int>& input)
{
  if(last-first > 1)
  {
    int mid = first + (last-first)/2;
    // cout << first << ", " << mid << ", " << last << "\n";
    #pragma omp parallel sections
    {
      #pragma omp section
      fragment_merge(first, mid, chunks, input);
      #pragma omp section
      fragment_merge(mid, last, chunks, input);
    }
    std::inplace_merge(input.begin() + chunks[first], input.begin() + chunks[mid], input.begin() + chunks[last]);
  }
}

int main(){
  
  FASTAFile file;
  vector<FASTAFile::FASTAEntry> entries;
  file.load("/buffer/ag_bsc/pmsb_23/data/sage/uniprot-SPTR_human_plusCont_FWBW.fasta", entries);
  //file.load("/buffer/ag_bsc/pmsb_23/max_alcer/sage/tests/Q99536.fasta", entries);
  cout << "Entires: " << entries.size() << "\n";
  
  double start;
  double end;

  start = omp_get_wtime();
  SearchDatabase sdb(entries);
  end = omp_get_wtime();

  cout << "construction time: " << end - start << "\n"; 

  MzMLFile f;

  MSExperiment exp;

  vector<pair<vector<SearchDatabase::Candidate>, size_t>> candidates;

  //vector<SearchDatabase::Candidate> candidates;

  f.load("/buffer/ag_bsc/pmsb_23/data/sage/centroided_Ms1andMs2.mzML", exp);
  //f.load("/buffer/ag_bsc/pmsb_23/max_alcer/sage/tests/LQSRPAAPPAPGPGQLTLR.mzML", exp);

  start = omp_get_wtime();
  sdb.search(exp, candidates);
  end = omp_get_wtime();

  cout << "search time: " << end - start << "\n";
  int sum = 0;
  for (auto& i : candidates){

    sum += i.first.size();

  }
  cout << "Erg: " << sum << "\n"; 
  /*

  SearchDatabase sdb({{"test", "test",
   "GSMTVDMQEIGSTEMPYEVPTQPNATSASAGRGWFDGPSFKVPSVPTRPSGIFRRPSRIKPEFSFKEKVSELVSPAVYTFGLFVQNASESLTSDDP"}});

  for (auto i : sdb.all_peptides_)
  {
    cout << i.sequence_ << "\n";
  }

  cout << "\n";

  for (auto h : sdb.bucket_frags_mz_)
  {
    cout << h << " ";
  }

  cout << "\n";
  cout << "\n";

  int counter = 0;

  for (auto j : sdb.all_fragments_)
  { 
    if(counter == sdb.bucketsize_)
    {
      cout << "\n";
      counter = 0;
    }

    cout << j.fragment_mz_ << " " << sdb.all_peptides_[j.peptide_index_].peptide_mz_ << "\n";

    counter++;
  }

  MSSpectrum spec;

  Precursor prec{};

  prec.setCharge(1);

  prec.setMZ(1039.48);

  spec.setPrecursors({prec});

  spec.push_back({391.176, 100});

  vector<SearchDatabase::Candidate> candidates;

  sdb.search(spec, candidates);

  for(auto c : candidates)
  {
    cout << *c.sequence_ << "\n ";
  }  
  */
  return 0;

}

