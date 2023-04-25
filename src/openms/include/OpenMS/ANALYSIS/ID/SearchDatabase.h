#pragma once

// #include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>

#include <vector>

namespace OpenMS
{

class OPENMS_DLLAPI SearchDatabase
{
    public:
    
    struct Candidate
    {   
        const AASequence* sequence_;
        size_t protein_index_; /// Index to Entry of std::vector<FASTAFile::FASTAEntry>
        Candidate(const AASequence* seq, size_t pi) : sequence_(seq), protein_index_(pi){}
    };

    SearchDatabase() = delete;
    /// Sets Default parameters and Builds up the Search Datastructure
    SearchDatabase(const std::vector<FASTAFile::FASTAEntry>& entries);  
    /// Searches Peaks of MSSpectrum in Data Structure 
    void search(MSExperiment& experiment, std::vector<std::pair<std::vector<Candidate>, size_t>>& candidates) const;

    void search(MSSpectrum& spectrum, std::vector<Candidate>& candidates) const;

    // private:

    // void updateMembers_() override;

    struct Fragment
    {
        size_t peptide_index_;
        double fragment_mz_;
        Fragment() = delete;
        Fragment(size_t prec, const Peak1D& frag):peptide_index_(prec), fragment_mz_(frag.getMZ()){}
    };

    struct Peptide
    {
        AASequence sequence_;
        size_t protein_index_;
        double peptide_mz_;
        Peptide() = default;
        Peptide(const Peptide&) = default;
        Peptide(const AASequence& pep, size_t index, double mass):sequence_(pep), protein_index_(index), peptide_mz_(mass){}
        Peptide& operator=(Peptide&&) = default;
        Peptide& operator=(const Peptide&) = default;
    };

    /// generates theoretical Peptides from all Proteins in fasta-File 
    std::vector<SearchDatabase::Peptide> generate_peptides_(const std::vector<FASTAFile::FASTAEntry>& entries) const;
    /// Merges presorted Chunks of Peptide-Fragments inplace
    void fragment_merge_(int first, int last, const std::vector<int>& chunks, std::vector<SearchDatabase::Fragment>& input) const;
    ///generates sortet vector with all theoretical Fragments for all theoretical Peptides
    std::vector<SearchDatabase::Fragment> generate_fragments_() const;

    std::string digestor_enzyme_;
    size_t missed_cleavages_;
    double peptide_min_mass_;
    double peptide_max_mass_;
    size_t peptide_min_length_;
    size_t peptide_max_length_;
    double fragment_min_mz_;
    double fragment_max_mz_;
    size_t bucketsize_;
    double precursor_mz_tolerance_;
    std::string precursor_mz_tolerance_unit_;
    double fragment_mz_tolerance_;
    std::string fragment_mz_tolerance_unit_;
    std::vector<Peptide> all_peptides_{};
    std::vector<double> bucket_frags_mz_{};
    std::vector<Fragment> all_fragments_{};

};

}