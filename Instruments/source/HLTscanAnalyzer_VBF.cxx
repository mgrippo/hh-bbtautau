/*! Merge DYJets files.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */
#include <boost/format.hpp>
#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "h-tautau/Analysis/include/TupleObjects.h"
#include "h-tautau/Analysis/include/SummaryTuple.h"
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
    TH2D_ENTRY(hlt_tau1_vs_hlt_tau2_2, 20, 20, 40, 20, 20, 40)
    TH2D_ENTRY(tau1_vs_tau2, 13, 27, 40, 13, 27, 40)
    TH1D_ENTRY(q2, 4, -2, 2)
    TH2D_ENTRY(hlt_tau1_vs_hlt_tau2_efficiency, 13, 27, 40, 13, 27, 40)
    TH2D_ENTRY(hlt_tau1_vs_hlt_tau2_efficiency_2, 20, 20, 40, 20, 20, 40)
        
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
    
    bool ApplyTauIdCut(const std::vector<uint32_t>& tauId_keys, const std::vector<float>& tauId_values) const
    {
        std::map<std::string, double> tau_id_cut;
        tau_id_cut.emplace("byMediumIsolationMVArun2v1DBoldDMwLT", 0.5);
        std::map<uint32_t, double> tau_id_cut_hashes;
        tau_id_cut_hashes.clear();
        for(const auto& id_cut: tau_id_cut)
            tau_id_cut_hashes[tools::hash(id_cut.first)] = id_cut.second;
        
        for(size_t n = 0; n < tauId_keys.size(); ++n) {
            if(tau_id_cut_hashes.count(tauId_keys.at(n)) &&
               tauId_values.at(n) < tau_id_cut_hashes.at(tauId_keys.at(n)))
                return false;
        }
        return true;
    }

    bool PassThreshold(const std::vector<std::pair<LorentzVectorE_Float, LorentzVectorE_Float>>& objects, double thr_hlt, double thr_l1)
    {
        bool pass = false;
        for(unsigned n = 0; n < objects.size(); ++n){
            std::pair<LorentzVectorE_Float, LorentzVectorE_Float> match = objects.at(n);
            LorentzVectorE_Float hlt_p4 = match.first;
            LorentzVectorE_Float l1_p4 = match.second;
            if( hlt_p4.Pt() >= thr_hlt && l1_p4.Pt() >= thr_l1)
                pass = true;
        }
        return pass;
    }
    
    void LoadInputs()
    {
        
        std::string filename = args.input_file();
        auto inputFile = root_ext::OpenRootFile(filename);
        ntuple::EventTuple eventTuple(args.tree_name(), inputFile.get(), true, {},
                                              GetEnabledBranches());
        const Long64_t n_entries = eventTuple.GetEntries();
        
        std::cout << "Entries: " << n_entries << std::endl;
        unsigned tot_events = 0;
        
        for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) { //loop on entries
            eventTuple.GetEntry(current_entry);
            const ntuple::Event& event = eventTuple.data();
            if (static_cast<EventEnergyScale>(event.eventEnergyScale) != analysis::EventEnergyScale::Central)
                continue;
            
            anaData.q2().Fill(event.q_2);
            LorentzVectorM_Float cand_tau_1;
            LorentzVectorM_Float cand_tau_2;
            
//            analysis::EventInfo<TauCandidate, TauCandidate> event_info(event, ntuple::JetPair{0, 1},nullptr);
//            
//            const auto& tau_1 = event_info.GetFirstLeg();
//            const auto& tau_2 = event_info.GetSecondLeg();
            
            if(!ApplyTauIdCut(event.tauId_keys_1, event.tauId_values_1) ||
               !ApplyTauIdCut(event.tauId_keys_2, event.tauId_values_2)) continue;
            
            if (event.extraelec_veto || event.extramuon_veto) continue;
            
            ++tot_events;
            
            //pair(HLT,L1)
            std::vector<std::pair<LorentzVectorE_Float,LorentzVectorE_Float>> hlt_l1_match_pair_1, hlt_l1_match_pair_2;
            
            for (unsigned n = 0; n < event.hlt_match_p4_1.size(); ++n){
                for(unsigned l = 0; l < event.l1_match_p4_1.size(); ++l){
                    LorentzVectorE_Float first_tau = event.hlt_match_p4_1.at(n);
                    LorentzVectorE_Float first_tau_l1 = event.l1_match_p4_1.at(l);
                    Float_t first_tau_l1_iso = event.l1_hwIso_1.at(l);
                    if (!(ROOT::Math::VectorUtil::DeltaR(first_tau,first_tau_l1) < 0.5)) continue;
                    if(first_tau_l1_iso < 0.5) continue;
                    std::pair<LorentzVectorE_Float,LorentzVectorE_Float> matched_pair_1(first_tau,first_tau_l1);
                    hlt_l1_match_pair_1.push_back(matched_pair_1);
                }
            }
            
            for (unsigned n = 0; n < event.hlt_match_p4_2.size(); ++n){
                for(unsigned l = 0; l < event.l1_match_p4_2.size(); ++l){
                    LorentzVectorE_Float second_tau = event.hlt_match_p4_2.at(n);
                    LorentzVectorE_Float second_tau_l1 = event.l1_match_p4_2.at(l);
                    Float_t second_tau_l1_iso = event.l1_hwIso_2.at(l);
                    if (!(ROOT::Math::VectorUtil::DeltaR(second_tau,second_tau_l1) < 0.5)) continue;
                    if(second_tau_l1_iso < 0.5) continue;
                    std::pair<LorentzVectorE_Float,LorentzVectorE_Float> matched_pair_2(second_tau,second_tau_l1);
                    hlt_l1_match_pair_2.push_back(matched_pair_2);
                }
            }
            
            
            for (unsigned n = 0; n <= anaData.hlt_tau1_vs_hlt_tau2().GetNbinsX(); ++n){
                for (unsigned h = 0; h <= anaData.hlt_tau1_vs_hlt_tau2().GetNbinsY(); ++h){
                    double value_x = anaData.hlt_tau1_vs_hlt_tau2().GetXaxis()->GetBinLowEdge(n);
                    double value_y = anaData.hlt_tau1_vs_hlt_tau2().GetYaxis()->GetBinLowEdge(h);
                    if (value_x > value_y) continue;
                    if(PassThreshold(hlt_l1_match_pair_1, value_x, value_x-2) && PassThreshold(hlt_l1_match_pair_2, value_y, value_y -2)){
                        double bin_content = anaData.hlt_tau1_vs_hlt_tau2().GetBinContent(n,h);
                        anaData.hlt_tau1_vs_hlt_tau2().SetBinContent(n,h,bin_content+1);
                    }
                }
            }
            
            for (unsigned n = 0; n <= anaData.hlt_tau1_vs_hlt_tau2_2().GetNbinsX(); ++n){
                for (unsigned h = 0; h <= anaData.hlt_tau1_vs_hlt_tau2_2().GetNbinsY(); ++h){
                    double value_x = anaData.hlt_tau1_vs_hlt_tau2_2().GetXaxis()->GetBinLowEdge(n);
                    double value_y = anaData.hlt_tau1_vs_hlt_tau2_2().GetYaxis()->GetBinLowEdge(h);
                    if (value_x > value_y) continue;
                    if(PassThreshold(hlt_l1_match_pair_1, value_x, value_x-2) && PassThreshold(hlt_l1_match_pair_2, value_y, value_y -2)){
                        double bin_content = anaData.hlt_tau1_vs_hlt_tau2_2().GetBinContent(n,h);
                        anaData.hlt_tau1_vs_hlt_tau2_2().SetBinContent(n,h,bin_content+1);
                    }
                }
            }
            

        } //end loop on entries

        double total_hlt_integral = anaData.hlt_tau1_vs_hlt_tau2().Integral();
        std::cout << "Total Integral: " << total_hlt_integral << std::endl;
        std::cout << "Tot events: " << tot_events << std::endl;
        
        anaData.hlt_tau1_vs_hlt_tau2_efficiency().CopyContent(anaData.hlt_tau1_vs_hlt_tau2());
        anaData.hlt_tau1_vs_hlt_tau2_efficiency().Scale(1./tot_events);
        
        anaData.hlt_tau1_vs_hlt_tau2_efficiency_2().CopyContent(anaData.hlt_tau1_vs_hlt_tau2_2());
        anaData.hlt_tau1_vs_hlt_tau2_efficiency_2().Scale(1./tot_events);
        
        
    }

};

} //namespace analysis

PROGRAM_MAIN(analysis::HLTscanAnalyzer, Arguments)
