#Loops_MC module

getItControllerName( itControllerType::Type{ItNumItController} ) = "itNumController";

getAttrLstItController( itControllerType::Type{ItNumItController} ) = ["itNum"];

getValLstItController( itController::ItNumItController ) = [itController.itNum];

testItNotDone( itController::ItNumItController ) = itController.itRef[] <= itController.itNum;

testItDoSample( itController::ItNumItController ) = (mod( itController.itRef[], itController.itSampleStep ) == 0);

testItDoStartSample( itController::ItNumItController ) = itController.itRef[] <= itController.itNumStartSample;

getItNumLst( itController::ItNumItController ) = (itController.itNum, itController.itNumSample, itController.itNumStartSample);




function getItControllerName( itControllerType::Type{JointItController} ) # = "itNumController";
	return join( ["joint", getItControllerName.(itControllerType.parameters[1].parameters)...], "_" );
end

function getAttrLstItController( itControllerType::Type{JointItController} )
	return append!( getAttrLstItController.( itControllerType.parameters[1].parameters ) );
end

function getValLstItController( itController::JointItController )
	return append!( getValLstItController.( itController.itControllerTup ) );
end

function testItNotDone( itController::JointItController )
	return reduce(&, testItNotDone.(itController.itControllerTup));
end

function testItDoSample( itController::JointItController )
	return reduce(|, testItDoSample.(itController.itControllerTup));
end

function testItDoStartSample( itController::JointItController )
	return reduce(|, testItDoStartSample.(itController.itControllerTup));
end

function getItNumLst( itController::JointItController )
	return zeros(Int64, 3);
end
