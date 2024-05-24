#Loops_MC module

function getInitializerName( initType::Type{BinomialInitializer} )
	return "BinomialInit";
end

function getAttrValInitializer( initializer::BinomialInitializer; rndDigs = rndDigsLpsMC )
	attrLst = [getInitializerName( initializer )];
	valLst = roundKeepInt.( [initializer.prob]; digits = rndDigs );
	
	return attrLst, valLst;
end

function initializeBL( initializer::BinomialInitializer, BfieldLst, params::ParamsLoops )
	for dim = 1 : params.nDimB
		rand!( initializer.dist, BfieldLst[dim] );
	end
end


function initializeBL( initializer::ConstantInitializer, BfieldLst, params::ParamsLoops )
	for dim = 1 : params.nDimB
		BfieldLst[dim] .= initializer.initVal;
	end
end

function getInitializerName( initType::Type{ConstantInitializer} )
	return "ConstInit";
end

function getAttrValInitializer( initializer::ConstantInitializer; rndDigs = rndDigsLpsMC )
	attrLst = [getInitializerName( initializer )];
	valLst = roundKeepInt.( [initializer.initVal]; digits = rndDigs );
	
	return attrLst, valLst;
end
