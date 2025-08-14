#ifndef INCLUDE_CONSTRUCT_TEST_DATA_H_
#define INCLUDE_CONSTRUCT_TEST_DATA_H_

#include "sdbg_build.h"
#include "settings.h"

Settings construct_default_settings() {
    Settings settings;

    settings.ram = 0.5;
    settings.threads = 1;
    settings.graph_folder = "graph";
    settings.cycles_folder = "cycles";
    settings.output_file = "CRISPR_Arrays.txt";

    return settings;
}

SDBG* construct_tiny_graph() {
    string tiny_graph = "test_data/tiny_graph/graph";

    SDBG* sdbg = new SDBG();
    sdbg->LoadFromFile(tiny_graph.c_str());

    return sdbg;
}

SDBG* construct_small_graph() {

}

SDBG* construct_medium_graph() {

}

#endif
