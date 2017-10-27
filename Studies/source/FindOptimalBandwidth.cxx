/*! Compute Optimal Bandwidth for each sample of mass or for ranges of masses.
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <initializer_list>
#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "AnalysisTools/Core/include/StatEstimators.h"
#include "hh-bbtautau/Analysis/include/MvaVariables.h"
#include "AnalysisTools/Core/include/NumericPrimitives.h"
#include "AnalysisTools/Core/include/ProgressReporter.h"
#include "AnalysisTools/Run/include/MultiThread.h"
#include "hh-bbtautau/Studies/include/MvaMethods.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Cuts/include/hh_bbtautau_2016.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "hh-bbtautau/Analysis/include/MvaConfigurationReader.h"

struct Arguments { // list of all program arguments
    REQ_ARG(std::string, input_path);
    REQ_ARG(std::string, cfg_file);
    REQ_ARG(unsigned, number_threads);
    REQ_ARG(bool, range);
    REQ_ARG(int, which_range);
    OPT_ARG(Long64_t, number_events, 5000);
};

namespace analysis {
namespace mva_study{

using clock = std::chrono::system_clock;

class FindOptimalBandwidth {
public:
    using Event = ntuple::Event;
    using EventTuple = ntuple::EventTuple;

    std::vector<ChannelSpin> set{{"muTau",0},{"eTau",0}, {"tauTau",0},{"muTau",2},{"eTau",2}, {"tauTau",2},
                                 {"tauTau",SM_spin}, {"muTau",SM_spin},{"eTau",SM_spin},
                                 {"muTau",bkg_spin},{"eTau",bkg_spin}, {"tauTau",bkg_spin}};

    FindOptimalBandwidth(const Arguments& _args): args(_args), vars(1, 12345678,{}, {"channel", "mass", "spin"}),
        reporter(std::make_shared<TimeReporter>())
    {
        MvaSetupCollection setups;
        SampleEntryListCollection samples_list;

        ConfigReader configReader;
        MvaConfigReader setupReader(setups);
        configReader.AddEntryReader("SETUP", setupReader, true);
        SampleConfigReader sampleReader(samples_list);
        configReader.AddEntryReader("FILES", sampleReader, false);
        configReader.ReadConfig(args.cfg_file());

        samples = samples_list.at("Samples").files;
    }

    void LoadSkimmedData()
    {
        for (const auto& s: set){
            std::cout << s.channel << s.spin <<std::endl;
            for(const SampleEntry& entry:samples)
            {
                if ( entry.spin != s.spin) continue;
                auto input_file = root_ext::OpenRootFile(args.input_path()+"/"+entry.filename);
                auto tuple = ntuple::CreateEventTuple(s.channel, input_file.get(), true, ntuple::TreeState::Skimmed);
                Long64_t tot_entries = 0;
                for(const Event& event : *tuple) {
                    if(tot_entries >= args.number_events()) break;
                    if (entry.id == SampleType::Bkg_TTbar && event.file_desc_id>=2) continue;
                    if (entry.id == SampleType::Sgn_NonRes && event.file_desc_id!=0) continue;
                    vars.AddEvent(event, entry.id, entry.spin, s.channel, entry.weight);
                    tot_entries++;
                }
                std::cout << entry << " number of events: " << tot_entries << std::endl;
            }
            samples_mass[s] = vars.GetSampleVariables(s.channel, s.spin);
            TimeReport();
        }
    }

    SampleIdVarData LoadRangeData(const SampleIdVarData& samples_mass){
        SampleIdVarData samples_range;
        for (const auto& range : analysis::mva_study::ranges){
            if (!range.Contains(args.which_range())) continue;
            for (const auto& sample : samples_mass){
                if (!range.Contains(sample.first.mass)) continue;
                for (const auto& entry : sample.second){
                    for (const auto& value: entry.second){
                        samples_range[SampleId{SampleType::Sgn_Res, range.min()}][entry.first].push_back(value);
                    }
                }
            }
        }
        return samples_range;
    }

    void TimeReport(bool tot = false) const
    {
        reporter->TimeReport(tot);
    }

    void Run()
    {
        run::ThreadPull threads(args.number_threads());
        LoadSkimmedData();

        std::string file_name_prefix = "OptimalBandwidth";
        if(args.range()) {
            std::cout<<"RANGES"<<std::endl;
            for (const auto& s: set){
                samples_range[s] = LoadRangeData(samples_mass[s]);
            }
            file_name_prefix += "Range";
        } else {
            samples_range = samples_mass;
        }

        for (const auto& s: set){
            std::cout<<s.channel<< "  " << s.spin<<std::endl;
            for (const auto& sample: samples_range[s]){
                std::cout<<"----"<<ToString(sample.first)<<"----"<<" entries: "<<sample.second.at("pt_l1").size()<<std::endl;
                bandwidth[s][sample.first] = OptimalBandwidth(sample.second);
                std::stringstream ss;
                ss << std::fixed << std::setprecision(0) << s.spin;
                std::string spin = ss.str();
                std::ofstream ListOptimalBandwidth(file_name_prefix+ToString(sample.first)+"_"+s.channel+"_spin"+spin+".csv", std::ofstream::out);
                for(const auto& value: bandwidth[s][sample.first]){
                   for(const auto& var_name: value.first)
                   {
                       ListOptimalBandwidth << var_name << "," ;
                   }
                   ListOptimalBandwidth << value.second << std::endl;
                }
                TimeReport();
            }
        }

        TimeReport(true);
    }

private:
    Arguments args;
    SampleEntryCollection samples;
    MvaVariablesStudy vars;
    std::shared_ptr<TimeReporter> reporter;
    std::vector<Range<int>> massranges;
    std::map<ChannelSpin,SampleIdVarData> samples_mass, samples_range;
    std::map<ChannelSpin,SampleIdNameElement> bandwidth, bandwidth_range;
};
}
}

PROGRAM_MAIN(analysis::mva_study::FindOptimalBandwidth, Arguments) // definition of the main program function