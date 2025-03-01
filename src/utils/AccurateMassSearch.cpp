// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTabMFile.h>
#include <OpenMS/FORMAT/OMSFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_AccurateMassSearch AccurateMassSearch

  @brief An algorithm to search for exact mass matches from a spectrum against a database (e.g. HMDB).

  <CENTER>
  <table>
  <tr>
  <th ALIGN = "center"> pot. predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=3> &rarr; AccurateMassSearch &rarr;</td>
  <th ALIGN = "center"> pot. successor tools </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderMetabo </td>
  <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> NA</td>
  </tr>
  </table>
  </CENTER>

  Accurate mass search against a database (usually HMDB).
  For details see @ref OpenMS::AccurateMassSearchEngine "AccurateMassSearchEngine".

  @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_AccurateMassSearch.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_AccurateMassSearch.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPAccurateMassSearch :
  public TOPPBase
{
public:
  TOPPAccurateMassSearch() :
    TOPPBase("AccurateMassSearch", "Match MS signals to molecules from a database by mass.", false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "featureXML or consensusXML file");
    setValidFormats_("in", {"featureXML", "consensusXML"});
    registerOutputFile_("out", "<file>", "", "mzTab file");
    setValidFormats_("out", ListUtils::create<String>("mzTab"));

    registerOutputFile_("out_annotation", "<file>", "", "A copy of the input file, annotated with matching hits from the database.", false);
    setValidFormats_("out_annotation", {"featureXML", "consensusXML", "oms"});

    // move some params from algorithm section to top level (to support input file functionality)
    Param p = AccurateMassSearchEngine().getDefaults();
    registerTOPPSubsection_("db", "Database files which contain the identifications");
    registerInputFileList_("db:mapping", "<file(s)>", ListUtils::toStringList<std::string>(p.getValue("db:mapping")), p.getDescription("db:mapping"), true, false, {"skipexists"});
    setValidFormats_("db:mapping", {"tsv"});
    registerInputFileList_("db:struct", "<file(s)>", ListUtils::toStringList<std::string>(p.getValue("db:struct")), p.getDescription("db:struct"), true, false, {"skipexists"});
    setValidFormats_("db:struct", {"tsv"});
    registerInputFile_("positive_adducts", "<file>", p.getValue("positive_adducts").toString(), p.getDescription("positive_adducts"), true, false, {"skipexists"});
    setValidFormats_("positive_adducts", {"tsv"});
    registerInputFile_("negative_adducts", "<file>", p.getValue("negative_adducts").toString(), p.getDescription("negative_adducts"), true, false, {"skipexists"});
    setValidFormats_("negative_adducts", {"tsv"});
    // addEmptyLine_();
    // addText_("Parameters for the accurate mass search can be given in the 'algorithm' part of INI file.");
    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    Param p = AccurateMassSearchEngine().getDefaults();
    // remove params which are already registered at top level (see registerOptionsAndFlags_())
    p.remove("db:mapping");
    p.remove("db:struct");
    p.remove("positive_adducts");
    p.remove("negative_adducts");
    return p;
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String file_ann = getStringOption_("out_annotation");

    Param ams_param = getParam_().copy("algorithm:", true);
    // copy top-level params to algorithm
    ams_param.setValue("db:mapping", ListUtils::create<std::string>(getStringList_("db:mapping")));
    ams_param.setValue("db:struct", ListUtils::create<std::string>(getStringList_("db:struct")));
    ams_param.setValue("positive_adducts", getStringOption_("positive_adducts"));
    ams_param.setValue("negative_adducts", getStringOption_("negative_adducts"));

    if (file_ann.hasSuffix("oms"))
    {
      ams_param.setValue("id_format", "ID"); // use IdentificationData to store id results
    }

    writeDebug_("Parameters passed to AccurateMassSearch", ams_param, 3);

    // mzTAB output data structure
    MzTab mztab_output;
    MzTabM mztabm_output;

    AccurateMassSearchEngine ams;
    ams.setParameters(ams_param);
    ams.init();

    std::string idf = std::string(ams.getParameters().getValue("id_format"));
    bool id_format = idf == "ID" ? true : false;

    FileTypes::Type filetype = FileHandler::getType(in);

    if (filetype == FileTypes::FEATUREXML)
    {
      FeatureMap ms_feat_map;
      FeatureXMLFile().load(in, ms_feat_map);

      //-------------------------------------------------------------
      // do the work
      //-------------------------------------------------------------
      if (id_format) // if format ID is used, MzTabM output will be generated
      {
        ams.run(ms_feat_map, mztabm_output);
      }
      else
      {
        ams.run(ms_feat_map, mztab_output);
      }

      //-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------

      if (file_ann.hasSuffix("featureXML"))
      {
        FeatureXMLFile().store(file_ann, ms_feat_map);
      }
      else if (file_ann.hasSuffix("oms"))
      {
        OMSFile().store(file_ann, ms_feat_map);
      }
    }
    else if (filetype == FileTypes::CONSENSUSXML && id_format)
    {
      throw Exception::InvalidValue(__FILE__,
                                    __LINE__,
                                    OPENMS_PRETTY_FUNCTION,
                                    "FATAL: CONSENSUSXML is currently not supporting ID and its MzTabM (v2.0.0-M) output, please use legacy_id",
                                    "");
    }
    else if (filetype == FileTypes::CONSENSUSXML)
    {
      ConsensusMap ms_cons_map;

      ConsensusXMLFile().load(in, ms_cons_map);

      //-------------------------------------------------------------
      // do the work
      //-------------------------------------------------------------
      ams.run(ms_cons_map, mztab_output);

      //-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------

      if (!file_ann.empty())
      {
        ConsensusXMLFile().store(file_ann, ms_cons_map);
      }
    }

    if(id_format && filetype == FileTypes::FEATUREXML)
    {
      MzTabMFile mztabm_file;
      mztabm_file.store(out, mztabm_output);
    }
    else
    {
      MzTabFile mztab_file;
      mztab_file.store(out, mztab_output);
    }

    return EXECUTION_OK;
  }
};


int main(int argc, const char** argv)
{
    TOPPAccurateMassSearch tool;
    return tool.main(argc, argv);
}

/// @endcond
