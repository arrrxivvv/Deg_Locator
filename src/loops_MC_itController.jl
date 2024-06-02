#Loops_MC module

getItControllerName( itControllerType::Type{ItNumItController} ) = "itNumController";

getAttrLstItController( itControllerType::Type{ItNumItController} ) = ["itNum"];

getValLstItController( itController::ItNumItController ) = [itController.itNum];

testItNotDone( itController::ItNumItController ) = itController.itRef[] <= itController.itNum;

testItDoSample( itController::ItNumItController ) = (mod( itController.itRef[], itController.itSampleStep ) == 0);

testItDoStartSample( itController::ItNumItController ) = itController.itRef[] <= itController.itNumStartSample;

getItNumLst( itController::ItNumItController ) = (itController.itNum, itController.itNumSample, itController.itNumStartSample);

