#include <OpenMS/ANALYSIS/ID/SearchDatabase.h>

//#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CHEMISTRY/AAIndex.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>

#include <omp.h>
#include <unordered_set>
#include <vector>

using namespace std;

namespace OpenMS
{
  
std::vector<SearchDatabase::Peptide> SearchDatabase::generate_peptides_(const std::vector<FASTAFile::FASTAEntry>& entries) const
{
  vector<SearchDatabase::Peptide> all_peptides;
  #pragma omp parallel
  { 
    ProteaseDigestion digestor;
    digestor.setEnzyme(digestor_enzyme_);
    digestor.setMissedCleavages(missed_cleavages_);
    vector<SearchDatabase::Peptide> all_peptides_pvt;
    #pragma omp for nowait
    for (size_t i = 0; i < entries.size(); i++)
    {
      vector<AASequence> peptides;
      digestor.digest(AASequence::fromString(entries[i].sequence), peptides, peptide_min_length_, peptide_max_length_);
      
      for (const auto& pep : peptides)
      { 
        if(pep.toString().find('X') != string::npos) continue;

        double seq_mz = pep.getMonoWeight();

        if(seq_mz < peptide_min_mass_ || seq_mz > peptide_max_mass_) continue;

        all_peptides_pvt.emplace_back(pep, i, seq_mz);
        
      }      

    }
    #pragma omp critical (add_critical)
    all_peptides.insert(all_peptides.end(), all_peptides_pvt.begin(), all_peptides_pvt.end());
  }
    return all_peptides;
}

void SearchDatabase::fragment_merge_(int first, int last, const std::vector<int>& chunks, std::vector<SearchDatabase::Fragment>& input) const
{
  if(last - first > 1)
  {
    int mid = first + (last-first)/2;
    // cout << first << ", " << mid << ", " << last << "\n";
    #pragma omp parallel sections
    {
      #pragma omp section
      SearchDatabase::fragment_merge_(first, mid, chunks, input);
      #pragma omp section
      SearchDatabase::fragment_merge_(mid, last, chunks, input);
    }
    std::inplace_merge(input.begin() + chunks[first], input.begin() + chunks[mid], input.begin() + chunks[last], 
    [](const  SearchDatabase::Fragment& l, const  SearchDatabase::Fragment& r)->bool
    {return (l.fragment_mz_ < r.fragment_mz_);});
  }
}

std::vector<SearchDatabase::Fragment> SearchDatabase::generate_fragments_() const
{
  TheoreticalSpectrumGenerator tsg;
  PeakSpectrum b_y_ions;
  std::vector<Fragment> all_frags;
  std::vector<int> chunk_start = {0};    
    
  for(size_t i = 0; i < all_peptides_.size(); i++)
  { 
    tsg.getSpectrum(b_y_ions, all_peptides_[i].sequence_, 1, 1);      
    for(const auto& frag : b_y_ions)
    { 
      if (frag.getMZ() < fragment_min_mz_ || frag.getMZ() > fragment_max_mz_) continue;
      all_frags.emplace_back(i, frag);        
    }
    chunk_start.emplace_back(all_frags.size());
    b_y_ions.clear(true);
  }

  SearchDatabase::fragment_merge_(0, chunk_start.size()-1, chunk_start, all_frags);

  cout << is_sorted(all_frags.begin(), all_frags.end(), 
  [](const  SearchDatabase::Fragment& l, const  SearchDatabase::Fragment& r)->bool
  {return (l.fragment_mz_ < r.fragment_mz_);}) << "\n";
   
  return all_frags;
}

SearchDatabase::SearchDatabase(const std::vector<FASTAFile::FASTAEntry>& entries)
{ 
  digestor_enzyme_ = "Trypsin";
  missed_cleavages_ = 0;
  peptide_min_mass_ = 500;
  peptide_max_mass_ = 5000;
  peptide_min_length_ = 5;
  peptide_max_length_ = 50;
  fragment_min_mz_ = 150;
  fragment_max_mz_ = 2000;
  bucketsize_ = 500000;
  precursor_mz_tolerance_ = 1;
  precursor_mz_tolerance_unit_ = "Da";
  fragment_mz_tolerance_ = 0.05;
  fragment_mz_tolerance_unit_ = "Da";

  all_peptides_ = generate_peptides_(entries);

  all_fragments_ = generate_fragments_();

  for(size_t i = 0; i < all_fragments_.size(); i += bucketsize_)
  {
    bucket_frags_mz_.emplace_back(all_fragments_[i].fragment_mz_);
  }
  #pragma omp parallel
  {
    #pragma omp for
    for(size_t i = 0; i < all_fragments_.size(); i += bucketsize_)
    {
      if (i+bucketsize_ > all_fragments_.size())
      {
        sort(all_fragments_.begin()+i, all_fragments_.end(), 
        [&](const SearchDatabase::Fragment& l, const SearchDatabase::Fragment& r)->bool 
        {return (all_peptides_[l.peptide_index_].peptide_mz_ < all_peptides_[r.peptide_index_].peptide_mz_);});
      }
      else
      {
        sort(all_fragments_.begin()+i, all_fragments_.begin()+i+bucketsize_, 
        [&](const SearchDatabase::Fragment& l, const SearchDatabase::Fragment& r)->bool 
        {return (all_peptides_[l.peptide_index_].peptide_mz_ < all_peptides_[r.peptide_index_].peptide_mz_);});
      }
    }
  }
  cout << "Buckets: " << bucket_frags_mz_.size() << "\n";

}

void SearchDatabase::search(MSSpectrum& spectrum, std::vector<SearchDatabase::Candidate>& candidates) const
{
  candidates.clear();
  unordered_set<size_t> index_hash;

  std::vector<Precursor> precursor = spectrum.getPrecursors();

  if(precursor.size() != 1) return;

  double prec_mz = precursor[0].getUnchargedMass();

  spectrum.sortByIntensity(true);
  size_t count_found = 0;
  size_t count_rounds = 0;

  for(const auto& peak : spectrum)
  {       
    if (count_found == 0 && count_rounds == 5) return;

    auto itr_lower = lower_bound(bucket_frags_mz_.begin(), bucket_frags_mz_.end(), peak.getMZ() - fragment_mz_tolerance_);
    
    size_t index_lower = distance(bucket_frags_mz_.begin(), itr_lower);
    if(index_lower != 0) index_lower--;
    while (bucket_frags_mz_[index_lower] <= (peak.getMZ() + fragment_mz_tolerance_) && index_lower < bucket_frags_mz_.size())
    {
      auto itr_lower_inner = lower_bound(all_fragments_.begin()+(index_lower*bucketsize_), all_fragments_.begin()+(index_lower*bucketsize_) + bucketsize_ -1, prec_mz - precursor_mz_tolerance_,
      [&](const SearchDatabase::Fragment& r, double l)->bool 
      {return (all_peptides_[r.peptide_index_].peptide_mz_ < l);});
      if(itr_lower_inner != all_fragments_.begin()+(index_lower*bucketsize_))--itr_lower_inner;
      auto itr_upper_inner = all_fragments_.begin() + (index_lower*bucketsize_) + bucketsize_ - 1;
      while(all_peptides_[(*itr_lower_inner).peptide_index_].peptide_mz_ <= (prec_mz + precursor_mz_tolerance_) && itr_lower_inner != itr_upper_inner)
      { 
        if(peak.getMZ() < (*itr_lower_inner).fragment_mz_ - fragment_mz_tolerance_ || peak.getMZ() > (*itr_lower_inner).fragment_mz_ + fragment_mz_tolerance_) continue;
        index_hash.insert((*itr_lower_inner).peptide_index_);
        ++itr_lower_inner;
        count_found++;
      }
      index_lower++;
    }
    cout << "\n";
    count_rounds++;
  }
  
  for (size_t j : index_hash)
  {
    candidates.emplace_back(&all_peptides_[j].sequence_, all_peptides_[j].protein_index_);
  }

}

void SearchDatabase::search(MSExperiment& experiment, std::vector<std::pair<std::vector<Candidate>, size_t>>& candidates) const
{
  candidates.clear();
  
  #pragma omp parallel
  { 
    #pragma omp for 
    for (size_t i = 0; i < experiment.size(); i++)
    { 
      unordered_set<size_t> index_hash;
      vector<Candidate> temp_cand;
      vector<Precursor> precursor = experiment[i].getPrecursors();

      if(precursor.size() != 1) continue;

      double prec_mz = precursor[0].getUnchargedMass();
      experiment[i].sortByIntensity(true);      
      size_t count_found = 0;
      size_t count_rounds = 0;

      for(const auto& peak : experiment[i])
      {        
        if (count_found == 0 && count_rounds == 5) break;

        auto itr_lower = lower_bound(bucket_frags_mz_.begin(), bucket_frags_mz_.end(), peak.getMZ() - fragment_mz_tolerance_);
        size_t index_lower = distance(bucket_frags_mz_.begin(), itr_lower);
        if(index_lower != 0) index_lower--;
        while (bucket_frags_mz_[index_lower] <= (peak.getMZ() + fragment_mz_tolerance_) && index_lower < bucket_frags_mz_.size())
        {
          auto itr_lower_inner = lower_bound(all_fragments_.begin()+(index_lower*bucketsize_), all_fragments_.begin()+(index_lower*bucketsize_) + bucketsize_ -1, prec_mz - precursor_mz_tolerance_,
          [&](const SearchDatabase::Fragment& r, double l)->bool {return (all_peptides_[r.peptide_index_].peptide_mz_ < l);});
          if(itr_lower_inner != all_fragments_.begin()+(index_lower*bucketsize_)) --itr_lower_inner;
          auto itr_upper_inner = all_fragments_.begin() + index_lower*bucketsize_ + bucketsize_-1;

          while(all_peptides_[(*itr_lower_inner).peptide_index_].peptide_mz_ <= (prec_mz + precursor_mz_tolerance_) && itr_lower_inner != itr_upper_inner)
          {   
            if((peak.getMZ() < (*itr_lower_inner).fragment_mz_ - fragment_mz_tolerance_) || (peak.getMZ() > (*itr_lower_inner).fragment_mz_ + fragment_mz_tolerance_))
            {
              ++itr_lower_inner;
              continue;
            }

            index_hash.insert((*itr_lower_inner).peptide_index_);
            ++itr_lower_inner;
            count_found++;
          }
          index_lower++;

        }
        count_rounds++;
      }      

      for (size_t j : index_hash)
      {
        temp_cand.emplace_back(&all_peptides_[j].sequence_, all_peptides_[j].protein_index_);
      }

      #pragma omp critical (add_critical)
      {
        candidates.emplace_back(temp_cand, i);
      }
    }
  }
}

} // ende namespace OpenMS





