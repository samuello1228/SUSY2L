#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ScanDir.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "SampleHandler/DiskListLocal.h"
#include <TSystem.h>
#include "SampleHandler/ScanDir.h"
#include "Analysis/TwoElectrons.h"

#include "EventLoop/LSFDriver.h"
#include "EventLoop/DirectDriver.h"
#include "EventLoopGrid/GridDriver.h"

#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoopAlgs/AlgSelect.h>

// Setup to take filename as argument?
int main (int argc, char* argv[] ) {
   
   bool batchMode = false;
   bool directMode = false;
   bool gridMode = false;
   
   bool limitEvents = false;
   int maxEvents = 0;
   int filesPerNode = 15;

   // Take the submit directory from the input if provided:
   std::string submitDir = "submit";
   if( argc > 1 ) submitDir = argv[ argc - 1 ];
   
   // Set up the job for xAOD access:
   xAOD::Init().ignore();

   // Construct the samples to run on:
   SH::SampleHandler sh;

   // use SampleHandler to open a list of files
   SH::readFileList(sh, "sample", "FileList.txt");

   
   // Set the name of the input TTree. It's always "CollectionTree"
   // for xAOD files.
   sh.setMetaString( "nc_tree", "CollectionTree" );

   // Print what we found:
   sh.print();

   // Create an EventLoop job:
   EL::Job job;
   job.sampleHandler( sh );

   // Add our analysis to the job:
   TwoElectrons* alg = new TwoElectrons;
   alg->submitDir = submitDir;
   
   
   // Handle command options
   if(argc > 2){
      for (int i = 1; i < argc - 1; i++){
         if(strcmp(argv[i], "-m") == 0) alg->matchMode = TwoElectrons::Matching::manualDR;
         else if (strcmp(argv[i], "-a") == 0) alg->matchMode = TwoElectrons::Matching::autoOld;
         else if (strcmp(argv[i], "-mc") == 0) alg->matchMode = TwoElectrons::Matching::MCTruthClass;
         
         else if(strcmp(argv[i], "-b") == 0) batchMode = true;
         else if(strcmp(argv[i], "-d") == 0) directMode = true;
         else if(strcmp(argv[i], "-g") == 0) gridMode = true;
         
         else if(strcmp(argv[i], "-f") == 0){
            i++;
            filesPerNode = atoi(argv[i]);
         }
         else if(strcmp(argv[i], "-n") == 0){
            i++;
            limitEvents = true;
            maxEvents = atoi(argv[i]);
         }
         else std::cout << "Invalid option: " << argv[i] << std::endl;
      }
   }
   // Default to directMode
   if((batchMode + directMode + gridMode) == 0){
      directMode = true;
   }
   
   // More than one selection
   if((batchMode + directMode + gridMode) != 1){
      std::cout << "ERROR:Please indicate -b batch mode, -d direct mode, or -g grid mode! \n\n";
      return 0;
   }
   
   if(alg->matchMode == TwoElectrons::Matching::manualDR){
      std::cout << "Manual MC truth matching.\n";
   } else if (alg->matchMode == TwoElectrons::Matching::autoOld){
      std::cout << "Auto (old) MC truth matching. \n";
   }
   else if (alg->matchMode == TwoElectrons::Matching::MCTruthClass){
      std::cout << "MCTruthClassifier matching. \n";
   }
   
   if(directMode) alg->testRunDirect = true;
   
   if(limitEvents){
      job.options()->setDouble (EL::Job::optMaxEvents, maxEvents);
   }
   
   job.algsAdd( alg );

   // Run the job on lxbatch
   if(batchMode){
      EL::LSFDriver driver;
      driver.options()->setString (EL::Job::optSubmitFlags, "-L /bin/bash"); // or whatever shell you are using
      driver.shellInit = "export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase && source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh";
      job.options()->setDouble (EL::Job::optFilesPerWorker, filesPerNode);
      
      std::cout << "Running on lxbatch with " << filesPerNode << " files per worker \n";
      
      driver.submit( job, submitDir );
   }
   if(directMode){
      EL::DirectDriver driver;
      
      std::cout << "Running directly on lxplus \n";
      
	   driver.submit( job, submitDir );
   }
   if(gridMode){
      EL::GridDriver driver;
      std::cout << "Submitting job to GRID" << std::endl;
      driver.submit( job, submitDir );
   }
   
   std::cout << std::endl;
   
   return 0;
   
}
