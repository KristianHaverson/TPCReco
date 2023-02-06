#include "EventSourceGRAW.h"
#include "EventSourceMultiGRAW.h"
#include "GeometryTPC.h"
#include <boost/program_options.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <iostream>

void eventProcessor(EventTPC *event) {
  auto *geometry = event->GetGeoPtr();
  // example
  // print waveform for cobo0, asad0, aget0, channel0
  // a shame we don't provide a standard iterator to access data without pain
  for (int i = 0; i != geometry->GetAgetNtimecells(); ++i) {
   // std::cout << event->GetValByAgetChannel(0, 0, 0, 0, i) << ",";
  }
  std::cout<<'\n';
}

void tpcDeco(const boost::property_tree::ptree &aConfig);

boost::program_options::variables_map parseCmdLineArgs(int argc, char **argv) {

  boost::program_options::options_description cmdLineOptDesc(
      "Allowed command line options");

  cmdLineOptDesc.add_options()("help", "produce help message")(
      "config,c", boost::program_options::value<std::string>()->required(),
      "string config file")(
      "dataFile,i", boost::program_options::value<std::string>(),
      "string - path to raw data file in single-GRAW mode (or list of "
      "comma-separated raw data files in multi-GRAW mode)")(
      "removePedestal", boost::program_options::value<bool>(),
      "bool - Flag to control pedestal removal. Overrides the value from "
      "config file.");

  boost::program_options::variables_map varMap;

  try {
    boost::program_options::store(
        boost::program_options::parse_command_line(argc, argv, cmdLineOptDesc),
        varMap);
    if (varMap.count("help")) {
      std::cout << std::endl
                << "tpcDeco config.json [options]" << std::endl
                << std::endl;
      std::cout << cmdLineOptDesc << std::endl;
      exit(1);
    }
    boost::program_options::notify(varMap);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    std::cout << cmdLineOptDesc << std::endl;
    exit(1);
  }

  return varMap;
}
/////////////////////////////////////
/////////////////////////////////////
int main(int argc, char **argv) {

  boost::program_options::variables_map varMap = parseCmdLineArgs(argc, argv);
  boost::property_tree::ptree tree;
  auto configFile = varMap["config"].as<std::string>();
  boost::property_tree::read_json(configFile, tree);

  if (varMap.count("dataFile")) {
    tree.put("dataFile", varMap["dataFile"].as<std::string>());
  }
  if (varMap.count("removePedestal")) {
    tree.put("removePedestal", varMap["removePedestal"].as<bool>());
  }
  if (tree.find("geometryFile") == tree.not_found()) {
    std::cerr << "missing geometryFile" << std::endl;
    return 1;
  }
  if (tree.find("dataFile") == tree.not_found()) {
    std::cerr << "missing geometryFile" << std::endl;
    return 1;
  }
  tpcDeco(tree);
  return 0;
}

void tpcDeco(const boost::property_tree::ptree &aConfig) {

  auto geometryFileName = aConfig.get<std::string>("geometryFile");
  auto dataFileName = aConfig.get<std::string>("dataFile");
  auto removePedestal = aConfig.get<bool>("removePedestal");
  auto singleAsadGrawFile = aConfig.get<bool>("singleAsadGrawFile");

  std::unique_ptr<EventSourceGRAW> myEventSource =
      singleAsadGrawFile
          ? std::make_unique<EventSourceMultiGRAW>(geometryFileName)
          : std::make_unique<EventSourceGRAW>(geometryFileName);
  myEventSource->setRemovePedestal(removePedestal);
  if (removePedestal) {
    if (aConfig.find("pedestal") != aConfig.not_found()) {
      myEventSource->configurePedestal(aConfig.find("pedestal")->second);
    }
  }
  myEventSource->loadDataFile(dataFileName);
  std::cout << "File with " << myEventSource->numberOfEntries()
            << " frames loaded." << std::endl;

  // initialize pedestal removal parameters for EventSource

  myEventSource->loadFileEntry(0);
  // finally acutal event loop
 // auto currentEventIdx = myEventSource->currentEventNumber();

/* do {

    auto event = myEventSource->getCurrentEvent();
    eventProcessor(event.get());
    // early return to skip reading full file in demo
      break;
    currentEventIdx = myEventSource->currentEventNumber();
    myEventSource->getNextEvent();
  } while (currentEventIdx != myEventSource->currentEventNumber());



*/

}
