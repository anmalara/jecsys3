#pragma once
#include "utils.hpp"

// static const std::vector<TString> samples = {"zjet", "gamjet", "hadw", "multijet", "incljet"};
static const std::vector<TString> samples = {"zjet", "gamjet", "hadw", "multijet"};
// static const std::vector<TString> samples = {"zjet", "gamjet", "hadw"};
// static const std::vector<TString> samples = {"zjet"};
// static const std::vector<TString> samples = {"gamjet"};
// static const std::vector<TString> samples = {"hadw"};
// static const std::vector<TString> samples = {"multijet"};
// static const std::vector<TString> samples = {"zjet", "gamjet"};
// static const std::vector<TString> samples = {"zjet", "multijet"};
// static const std::vector<TString> samples = {"zjet", "hadw"};
// static const std::vector<TString> samples = {"gamjet", "hadw"};
// static const std::vector<TString> samples = {"gamjet", "multijet"};
// static const std::vector<TString> samples = {"hadw", "multijet"};

static const std::vector<TString> hdm_methods = {"mpf", "db"};
// static const std::vector<TString> hdm_methods = {"mpf"};
// static const std::vector<TString> types = {"Resp","chf","nef","nhf"};
static const std::vector<TString> types = {"Resp"};

static const TString output_fname = "rootfiles/outputRUN.root";

static constexpr double ptmax_multijet = 300;
