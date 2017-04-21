/*! Merge DYJets files.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"


struct Arguments {
    run::Argument<std::string> tree_name{"tree_name", "Tree on which we work"};
    run::Argument<std::string> input_file{"input_file", "Input file of the sample with full path"};
    run::Argument<std::string> output_file{"output_file", "Output file"};
};

namespace analysis {
    
class HLTscanAnalyzerData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;
    
    //TH2D_ENTRY(lhe_hh_cosTheta_vs_m, 25, 200, 2000, 30, -1.1, 1.1)
        
};
    
    
class HLTscanAnalyzer {
public:
    HLTscanAnalyzer(const Arguments& _args) :
    args(_args), output(root_ext::CreateRootFile(args.output_file())), anaData(output)
    {
        LoadInputs();
    }
        
    void Run() {}
        
private:
    Arguments args;
    std::shared_ptr<TFile> output;
    HLTscanAnalyzerData anaData;


    static const std::set<std::string>& GetEnabledBranches()
    {
        static const std::set<std::string> EnabledBranches_read = {
            "l1_match_p4_1", "l1_match_p4_2", "l1_hwIso_1", "l1_hwIso_2",
            "hlt_match_p4_1", "hlt_match_p4_2"
        };
        return EnabledBranches_read;
    }

    
    void LoadInputs()
    {
        std::string filename = args.input_file();
        auto inputFile = root_ext::OpenRootFile(filename);
        ntuple::EventTuple eventTuple(args.tree_name(), inputFile.get(), true, {},
                                              GetEnabledBranches());
        const Long64_t n_entries = eventTuple.GetEntries();

        for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) { //loop on entries
            eventTuple.GetEntry(current_entry);
            const ntuple::Event& event = eventTuple.data();
            if (static_cast<EventEnergyScale>(event.eventEnergyScale) != analysis::EventEnergyScale::Central)
                continue;
            std::cout << "l1 hwIso 1:" << event.l1_hwIso_1 << std::endl;
            std::cout << "l1 hwIso 2:" << event.l1_hwIso_2 << std::endl;
        } //end loop on entries
    }

};

} //namespace analysis

PROGRAM_MAIN(analysis::HLTscanAnalyzer, Arguments)
