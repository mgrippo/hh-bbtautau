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
    TH2D_ENTRY(hlt_tau1_vs_hlt_tau2, 13, 27, 40, 13, 27, 40)
    TH2D_ENTRY(tau1_vs_tau2, 13, 27, 40, 13, 27, 40)
    TH1D_ENTRY(q2, 4, -2, 2)
    TH2D_ENTRY(hlt_tau1_vs_hlt_tau2_efficiency, 13, 27, 40, 13, 27, 40)
        
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
            "p4_1", "p4_2","q_2", "l1_match_p4_1", "l1_match_p4_2", "l1_hwIso_1", "l1_hwIso_2",
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
//            for (unsigned n = 0; n < event.l1_hwIso_1.size(); ++n){
//                std::cout << "l1 hwIso 1:" << event.l1_hwIso_1.at(n) << std::endl;
//            }
//            for (unsigned n = 0; n < event.l1_hwIso_2.size(); ++n){
//                std::cout << "l1 hwIso 2:" << event.l1_hwIso_2.at(n) << std::endl;
//            }

            anaData.q2().Fill(event.q_2);
            for (unsigned n = 0; n < event.hlt_match_p4_1.size(); ++n){
                for (unsigned h = 0; h < event.hlt_match_p4_2.size(); ++h){
                    LorentzVectorE_Float first_tau = event.hlt_match_p4_1.at(n);
                    LorentzVectorE_Float second_tau = event.hlt_match_p4_2.at(h);
                    if (first_tau.Pt() >= 27 && first_tau.Pt() <= 40){
                        if (second_tau.Pt() >= 27 && second_tau.Pt() <= 40){
                            anaData.hlt_tau1_vs_hlt_tau2().Fill(second_tau.Pt(), first_tau.Pt());
                        }
                    }
                }
            }


            if (event.p4_1.pt() >= 27 && event.p4_1.pt() <= 40){
                if (event.p4_2.pt() >= 27 && event.p4_2.pt() <= 40){
                    anaData.tau1_vs_tau2().Fill(event.p4_2.pt(), event.p4_1.pt());
                }
            }

        } //end loop on entries

        double total_hlt_integral = anaData.hlt_tau1_vs_hlt_tau2().Integral();
        std::cout << "Total Integral: " << total_hlt_integral << std::endl;
        double total_integral = anaData.tau1_vs_tau2().Integral();
        std::cout << "Total Integral: " << total_integral << std::endl;
        for (unsigned n = 0; n < anaData.hlt_tau1_vs_hlt_tau2().GetNbinsX(); ++n){
            for (unsigned h = 0; h < anaData.hlt_tau1_vs_hlt_tau2().GetNbinsY(); ++h){
                double bin_content = anaData.hlt_tau1_vs_hlt_tau2().GetBinContent(n,h);
                double bin_efficiency = bin_content/total_integral;
                anaData.hlt_tau1_vs_hlt_tau2_efficiency().SetBinContent (n,h,bin_efficiency);
            }
        }
    }

};

} //namespace analysis

PROGRAM_MAIN(analysis::HLTscanAnalyzer, Arguments)
