#include "InputFileHelper.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <boost/filesystem.hpp>
#include <fstream>
#include <set>

namespace fs = boost::filesystem;
using namespace InputFileHelper;
class InputFileHelperTest : public ::testing::Test {
public:
  static std::string directory;

  static std::string createFile(const std::string &name) {
    auto file = directory + name;
    std::ofstream{file};
    return file;
  }

  static void TearDownTestSuite() { fs::remove_all(directory); }

  static void SetUpTestSuite() {
    directory = (fs::temp_directory_path() / fs::unique_path()).string() +
                fs::path::preferred_separator;
    fs::create_directories(directory);

    createFile("CoBo_ALL_AsAd_ALL_2021-07-12T11:02:15.328_0000.graw");
    createFile("CoBo_ALL_AsAd_ALL_2021-07-12T11:02:15.328_0001.graw");
    createFile("CoBo_ALL_AsAd_ALL_2021-07-12T11:02:15.328_0002.graw");

    createFile("CoBo0_AsAd0_2021-07-12T12:03:40.978_0000.graw");
    createFile("CoBo0_AsAd1_2021-07-12T12:03:40.982_0000.graw");
    createFile("CoBo0_AsAd2_2021-07-12T12:03:40.994_0000.graw");
    createFile("CoBo0_AsAd3_2021-07-12T12:03:41.001_0000.graw");

    createFile("CoBo0_AsAd0_2021-07-12T12:03:40.978_0001.graw");
    createFile("CoBo0_AsAd1_2021-07-12T12:03:40.982_0001.graw");
    createFile("CoBo0_AsAd2_2021-07-12T12:03:40.994_0001.graw");
    createFile("CoBo0_AsAd3_2021-07-12T12:03:41.001_0001.graw");

    createFile("1.txt");
    createFile("1.graw");
  }
};

std::string InputFileHelperTest::directory = "";

TEST_F(InputFileHelperTest, Tokenizer) {
  EXPECT_THAT(
      tokenize("CoBo0_AsAd0_2021-07-12T12:03:40.978_0000.graw,CoBo0_AsAd1_"
               "2021-07-12T12:03:40.982_0000.graw"),
      ::testing::ElementsAreArray(
          {"CoBo0_AsAd0_2021-07-12T12:03:40.978_0000.graw",
           "CoBo0_AsAd1_2021-07-12T12:03:40.982_0000.graw"}));
  std::vector<std::string> token = {
      "CoBo0_AsAd0_2021-07-12T12:03:40.978_0000.graw"};
  EXPECT_THAT(
      tokenize("CoBo0_AsAd0_2021-07-12T12:03:40.978_0000.graw"),
      ::testing::ElementsAre("CoBo0_AsAd0_2021-07-12T12:03:40.978_0000.graw"));
}

TEST_F(InputFileHelperTest, ExtensionValidation) {
  {
    std::vector<std::string> files = {
        "CoBo_ALL_AsAd_ALL_2021-07-12T11:02:15.328_0000.graw",
        "CoBo_ALL_AsAd_ALL_2021-07-12T11:02:15.328_0001."
        "root"};
    EXPECT_THROW(getExtension(files.begin(), files.end()), std::runtime_error);
  }
  {
    std::vector<std::string> files = {
        "CoBo_ALL_AsAd_ALL_2021-07-12T11:02:15.328_0000.graw",
        "CoBo_ALL_AsAd_ALL_2021-07-12T11:02:15.328_0001.graw"};
    EXPECT_EQ(getExtension(files.begin(), files.end()), ".graw");
  }
}

TEST_F(InputFileHelperTest, FileDiscovery) {
  std::vector<std::string> files;
  discoverFiles(directory + "CoBo0_AsAd0_2021-07-12T12:03:40.978_0000.graw",
                std::chrono::milliseconds(30), std::back_inserter(files));
  EXPECT_THAT(
      files,
      ::testing::UnorderedElementsAreArray(
          {directory + "CoBo0_AsAd0_2021-07-12T12:03:40.978_0000.graw",
           directory + "CoBo0_AsAd1_2021-07-12T12:03:40.982_0000.graw",
           directory + "CoBo0_AsAd2_2021-07-12T12:03:40.994_0000.graw",
           directory + "CoBo0_AsAd3_2021-07-12T12:03:41.001_0000.graw"}));
}

TEST_F(InputFileHelperTest, FileDiscoveryCSV) {
  auto files = discoverFilesCSV(
      directory + "CoBo0_AsAd0_2021-07-12T12:03:40.978_0000.graw",
      std::chrono::milliseconds(30));
  EXPECT_EQ(files,
            directory + "CoBo0_AsAd0_2021-07-12T12:03:40.978_0000.graw" + "," +
                directory + "CoBo0_AsAd1_2021-07-12T12:03:40.982_0000.graw" +
                "," + directory +
                "CoBo0_AsAd2_2021-07-12T12:03:40.994_0000.graw" + "," +
                directory + "CoBo0_AsAd3_2021-07-12T12:03:41.001_0000.graw");
}

TEST(InputFileHelper, ExtensionsFilter) {
  std::set<std::string> allowedExtensions = {".graw", ".root", ".log"};
  std::vector<std::string> files = {
      "CoBo0_AsAd0_2021-07-12T12:03:40.978_0000.graw",
      "CoBo0_AsAd0_2021-07-12T12:03:40.978_0000.root",
      "CoBo0_AsAd0_2021-07-12T12:03:40.978_0000.graw.log",
      "CoBo0_AsAd0_2021-07-12T12:03:40.978_0000",
      "CoBo0_AsAd0_2021-07-12T12:03:40.978_0000.txt",
      "CoBo0_AsAd0_2021-07-12T12:03:40.978_0000.graw.txt"};
  auto last =
      filterExtensions(std::begin(files), std::end(files), allowedExtensions);
  files.erase(last, std::end(files));
  EXPECT_THAT(files,
              ::testing::UnorderedElementsAreArray(
                  {"CoBo0_AsAd0_2021-07-12T12:03:40.978_0000.graw",
                   "CoBo0_AsAd0_2021-07-12T12:03:40.978_0000.root",
                   "CoBo0_AsAd0_2021-07-12T12:03:40.978_0000.graw.log"}));
}

TEST_F(InputFileHelperTest, FileDiscoveryFromRunId) {
  std::vector<std::string> files;
  std::string runId = "20210712120340";
  auto timePoint = RunIdParser::timePointFromRunId(runId);
  unsigned long chunk = 0;
  auto begin =
      boost::filesystem::directory_iterator(boost::filesystem::path(directory));
  auto end = boost::filesystem::directory_iterator{};

  discoverFiles(timePoint, chunk, std::chrono::milliseconds(1500), begin, end,
                std::back_inserter(files));
  EXPECT_THAT(
      files,
      ::testing::UnorderedElementsAreArray(
          {directory + "CoBo0_AsAd0_2021-07-12T12:03:40.978_0000.graw",
           directory + "CoBo0_AsAd1_2021-07-12T12:03:40.982_0000.graw",
           directory + "CoBo0_AsAd2_2021-07-12T12:03:40.994_0000.graw",
           directory + "CoBo0_AsAd3_2021-07-12T12:03:41.001_0000.graw"}));
  files.clear();
  begin =
      boost::filesystem::directory_iterator(boost::filesystem::path(directory));
  discoverFiles(timePoint, chunk, std::chrono::milliseconds(1000), begin, end,
                std::back_inserter(files));
  EXPECT_THAT(files,
              ::testing::UnorderedElementsAreArray({
                  directory + "CoBo0_AsAd0_2021-07-12T12:03:40.978_0000.graw",
                  directory + "CoBo0_AsAd1_2021-07-12T12:03:40.982_0000.graw",
                  directory + "CoBo0_AsAd2_2021-07-12T12:03:40.994_0000.graw"
                  // last file missing becasue not in duration range
              }));
}
