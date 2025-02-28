from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from MSSpectrum cimport *
from AASequence cimport *
from ProteinProteinCrossLink cimport *

cdef extern from "<OpenMS/CHEMISTRY/SimpleTSGXLMS.h>" namespace "OpenMS":

    cdef cppclass SimpleTSGXLMS(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler

        SimpleTSGXLMS() nogil except +
                # wrap-doc:
                #  Generates theoretical spectra for cross-linked peptides
                #  
                #  The spectra this class generates are vectors of SimplePeaks
                #  This class generates the same peak types as TheoreticalSpectrumGeneratorXLMS
                #  and the interface is very similar, but it is simpler and faster
                #  SimplePeak only contains an mz value and a charge. No intensity values
                #  or String annotations or other additional DataArrays are generated

        SimpleTSGXLMS(SimpleTSGXLMS &) nogil except +

        void getLinearIonSpectrum(libcpp_vector[ SimplePeak ]& spectrum, AASequence peptide,
                Size link_pos, int charge, Size link_pos_2) nogil except +
                # wrap-doc:
                #  Generates fragment ions not containing the cross-linker for one peptide
                #  
                #  B-ions are generated from the beginning of the peptide up to the first linked position,
                #  y-ions are generated from the second linked position up the end of the peptide
                #  If link_pos_2 is 0, a mono-link or cross-link is assumed and the second position is the same as the first position
                #  For a loop-link two different positions can be set and link_pos_2 must be larger than link_pos
                #  The generated ion types and other additional settings are determined by the tool parameters
                #  
                #  :param spectrum: The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it
                #  :param peptide: The peptide to fragment
                #  :param link_pos: The position of the cross-linker on the given peptide
                #  :param charge: The maximal charge of the ions
                #  :param link_pos_2: A second position for the linker, in case it is a loop link

        void getXLinkIonSpectrum(libcpp_vector[ SimplePeak ]& spectrum, AASequence peptide,
                Size link_pos, double precursor_mass,
                int mincharge, int maxcharge, Size link_pos_2) nogil except +
                # wrap-doc:
                #  Generates fragment ions containing the cross-linker for one peptide
                #  
                #  B-ions are generated from the first linked position up to the end of the peptide,
                #  y-ions are generated from the beginning of the peptide up to the second linked position
                #  If link_pos_2 is 0, a mono-link or cross-link is assumed and the second position is the same as the first position
                #  For a loop-link two different positions can be set and link_pos_2 must be larger than link_pos
                #  Since in the case of a cross-link a whole second peptide is attached to the other side of the cross-link,
                #  a precursor mass for the two peptides and the linker is needed
                #  In the case of a loop link the precursor mass is the mass of the only peptide and the linker
                #  Although this function is more general, currently it is mainly used for loop-links and mono-links,
                #  because residues in the second, unknown peptide cannot be considered for possible neutral losses
                #  The generated ion types and other additional settings are determined by the tool parameters
                #  
		#  :param spectrum: The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it
		#  :param peptide: The peptide to fragment
		#  :param link_pos: The position of the cross-linker on the given peptide
		#  :param precursor_mass: The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum
		#  :param mincharge: The minimal charge of the ions
		#  :param maxcharge: The maximal charge of the ions, it should be the precursor charge and is used to generate precursor ion peaks
		#  :param link_pos_2: A second position for the linker, in case it is a loop link

        void getXLinkIonSpectrum(libcpp_vector[ SimplePeak ]& spectrum, ProteinProteinCrossLink crosslink,
                bool frag_alpha, int mincharge, int maxcharge) nogil except +
                # wrap-doc:
                #  Generates fragment ions containing the cross-linker for a pair of peptides
                #  
                #  B-ions are generated from the first linked position up to the end of the peptide,
                #  y-ions are generated from the beginning of the peptide up to the second linked position
                #  This function generates neutral loss ions by considering both linked peptides
                #  Only one of the peptides, decided by @frag_alpha, is fragmented
                #  This simplifies the function, but it has to be called twice to get all fragments of a peptide pair
                #  The generated ion types and other additional settings are determined by the tool parameters
                #  This function is not suitable to generate fragments for mono-links or loop-links
                #  
		#  :param spectrum: The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it
                #  :param crosslink: ProteinProteinCrossLink to be fragmented
                #  :param link_pos: The position of the cross-linker on the given peptide
                #  :param precursor_mass: The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum
                #  :param frag_alpha: True, if the fragmented peptide is the Alpha peptide
                #  :param mincharge: The minimal charge of the ions
                #  :param maxcharge: The maximal charge of the ions, it should be the precursor charge and is used to generate precursor ion peaks


cdef extern from "<OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>" namespace "OpenMS::SimpleTSGXLMS":

    cdef cppclass SimplePeak "OpenMS::SimpleTSGXLMS::SimplePeak":

        SimplePeak() nogil except + # wrap-doc:A simple struct to represent peaks with mz and charge and sort them easily
        SimplePeak(double mz, int charge) nogil except +
        SimplePeak(SimplePeak &) nogil except +

        double mz
        int charge
