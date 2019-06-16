#!/bin/bash
fake_subMit_command fake_computeNodes fake_memory fakeTime fake_jobName fakeOut.o
mkdir -p fake_serverPath/fakeDir/
cd fake_workDir/ ; cp *.sqlite *.RData fake_serverPath/fakeDir/
fake_appPath fake_commandArgs fake_passedArgs fake_workDir/fake_tmpScript.r fake_workDir/fake_tmpScript.r.Rout
cd fake_serverPath/fakeDir/ ; mv *.sqlite *.Rout *.RData fake_workDir/
