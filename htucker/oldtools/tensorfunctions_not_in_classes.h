#ifndef LAWA_HTUCKER_TENSORFUNCTIONS_H
#define LAWA_HTUCKER_TENSORFUNCTIONS_H 1

#include <fima/HTucker/dimensionindex.h>

double inverseNorm(DimensionIndex para){
	double res = 0.0;
	for(int i = 0; i < para.length(); ++i){
		res+= para[i]*para[i];
	}
	return 1.0/sqrt(res);
};
  
double binary(DimensionIndex para){
	double res = 1.0;
	for(int i = 0; i < para.length(); ++i){
		res+= pow(2.0,i)*(para[i]-1);
	}
	return res;
};

double decimal(DimensionIndex para){
	double res = 1.0;
	for(int i = 0; i < para.length(); ++i){
		res+= pow(10.0,i)*(para[i]-1);
	}
	return res;
};

	
void  getVecsexpsumApproximation(DenseVectorList<double> & list,int d,DimensionIndex minvec, DimensionIndex maxvec){
	DenseVector<Array<double> > weights(40);
	weights =  0.0004541105141216453968248119678948609979,
			   0.0078530680104171971521949987930788505963,
			 0.0358442951018547386594730302239453223478,
			 0.0800511636731031012202772917385473760987,
			 0.1168715918274595083080672451580905857327,
			 0.1317468847713531808620935140141661179314 ,
			 0.1268550732122552672878103285913908848670 ,
			 0.1109657824442714362255012552904709366430,
			 0.0915407675164530429816095620509930341768,
			 0.0728296028697006614068595181654686143702,
			 0.0566262616355938613997447837050680519155,
			 0.0433600894479567821162072083218674833915,
			 0.0328421513638694122585743456754769908912,
			 0.0246657172704435717145978130063199706967,
			 0.0183917192322462027550102406370235996746 ,
			 0.0136227224770988644418135728175744958435,
			 0.0100250532523774125004520276632780095838,
			 0.0073290867029904425973310905066682163778,
			 0.0053216028826095790999957683888726966970 ,
			 0.0038362256666861645005019393166179053622,
			 0.0027443224789343381565523340582385403152,
			 0.0019471510649087046089954599194064843726,
			 0.0013693829522368219844082718491921901727,
			 0.0009538877370448752111114987457026609396,
			 0.0006575902435281989381065818672453260341,
			 0.0004482091681493420014885365430801211772 ,
			 0.0003017066922888354188646730291625125409,
			 0.0002003052463432384262700059745642040820,
			 0.0001309533404038651058759093782370673542,
			 0.0000841450289679689093394985195683861284,
			 0.0000530163906346343672674381358229543301,
			 0.0000326579433267181342724703319572348084,
			 0.0000195946088007008467557055835143813537,
			 0.0000113949728253805224946867206361372001,
			 0.0000063797866387841876792628796961806114,
			 0.0000034061919433246302751298561507290710,
			 0.0000017093587076932593639466942818314337,
			 0.0000007873796897325934882231604846542004,
			 0.0000003186056947254193096300954077315363,
			 0.0000001024294180134919968361669149146366;
	DenseVector<Array<double> > alpha(40);
	alpha = 0.0328758165954420368502493336737568796480,
			 0.0570473982762439551973880462465427143570,
			 0.0966283152617544954693135257628444634292,
			 0.1636063258412157908437819728297846211262,
			 0.2786723770321616242787599732055880963344,
			 0.4784495396399847536917683543666868217770,
			 0.8285231470543625377543975329963643616793,
			 1.4475289223590702735247615073888027836801,
			 2.5522100303561091739872057626570267530042,
			 4.5426301651741418010516104786233881895896,
			 8.1651768488703318207497239278325196210062 ,
			14.8281159280173661825508180100996469263919,
			27.2202001694938288487574107321620431321207,
			50.5395968142454580715461354145645600510761,
			94.9698324334290623002785558526284148683771,
			180.7430429547675229534631213823558937292546,
			348.6595709978226785563659717581685981713235,
			682.3160325424825219386981700608885148540139 ,
			1355.9285760062651674529021761372860055416822,
			2739.2321520142396620300928589131217449903488,
			5632.4067005890412889179685862472979351878166,
			11804.1327408178242155400994306546635925769806,
			25254.1956632226066332691516436170786619186401,
			55255.5977602627446430005875299684703350067139,
			123898.7837399924441115217632614076137542724609,
			285403.3681506335250048778107156977057456970215 ,
			677314.9416636431666347561986185610294342041016,
			1661627.5164404432468927552690729498863220214844,
			4231169.6896806240947626065462827682495117187500,
			11239234.2676030291650022263638675212860107421875,
			31336866.8801911466143792495131492614746093750000,
			92434878.3605544002202805131673812866210937500000,
			291426049.6279058846412226557731628417968750000000,
			995657663.5402839354937896132469177246093750000000 ,
			3757776658.0723961133044213056564331054687500000000,
			16116802909.0660596201196312904357910156250000000000,
			82144510740.0004426687955856323242187500000000000000,
			538206356175.0550399124622344970703125000000000000000,
			5341579170601.2921476364135742187500000000000000000000,
			130450953442068.3405914306640625000000000000000000000000;

	
				for(int i = 0; i< d; i++){
					for(int j = 1; j<= 40; j++){
						DenseVector<Array<double> > tmp(maxvec[i]-minvec[i]+1);
						for(int l= tmp.firstIndex(); l<=tmp.lastIndex(); ++l){
							//std::cout << "i= " << i << "  j= " << j << "    " << (minvec(i)+l-tmp.firstIndex()) << endl;
							if(i == 1){
								tmp(l) = weights(j)*exp(-((minvec[i]+l-tmp.firstIndex())/(maxvec[i]-minvec[i]+1.0))*((minvec[i]+l-tmp.firstIndex())/(maxvec[i]-minvec[i]+1.0))*alpha(j));
							} else {
								tmp(l) = exp(-((minvec[i]+l-tmp.firstIndex() + 0.0)*(minvec[i]+l-tmp.firstIndex()))/((maxvec[i]-minvec[i]+1.0)*(maxvec[i]-minvec[i]+1.0) + 0.0)*alpha(j));
							
							}

						}
						list.add(tmp);
					}
				}


};

void getVecsinverseNormApproximation(DenseVectorList<double> & list,int d,DimensionIndex minvec, DimensionIndex maxvec){
	assert(d == minvec.length() && d == maxvec.length());
	DenseVector<Array<double> > weights(35);
	weights =  0.0003476617162031463656887349288667626857,
		 0.0003656677800731512040516686207159910427,
		 0.0004075639038080694722811763840587123697,
		 0.0004843102689961686493529664332240158853,
		 0.0006075484362407267528193426119921538575,
		 0.0007860621238313433964849066446726327806,
		 0.0010286920148215165835244833104043793437,
		 0.0013484039992954368087911095805929428959,
		 0.0017637932498565288986423389163381292288,
		 0.0022995948323400161628149234675083656398,
		 0.0029874656853097934646253271400578865880,
		 0.0038672614791946895807455515312045513099,
		 0.0049887848083617705673121177356588162866,
		 0.0064140049762051946081513006427121670328,
		 0.0082197972649914829865058925849330417890,
		 0.0105012852751240196090380674132092053696,
		 0.0133758979709688672253838190480401681981,
		 0.0169882797026389577249398925989920527968,
		 0.0215162197131539093567588515412691840822,
		 0.0271778000398019895103514419545942537582,
		 0.0342399988121945339294440423250076577233,
		 0.0430290321349169279780387271233665913428,
		 0.0539427787599617519049503139294543530013,
		 0.0674657211055092835482404310665227598065,
		 0.0841869982556718884943842105827993549383,
		 0.1048225290121168909272938649102080432840,
		 0.1302430953236043887537645713514677936473,
		 0.1615128713218071510428555659366622876405,
		 0.1999503386727134369738609689926001067306,
		 0.2472453213059679058203037954610792326093,
		 0.3057312880122441022998601434235155238639,
		 0.3791210128165055750809782353361043760742,
		 0.4747822440254026090267897020069653990504,
		 0.6124014951064692475124137671649293679366,
		 0.8770631816903148868821808725737554368607;
	DenseVector<Array<double> > alpha(35);
	alpha = 0.0000000235400856956305055878655722424688,
		 0.0000002190687866559873968066208739485593 ,
		 0.0000006537664626626933778312767948510514,
		 0.0000014417390459603761640437794711739619,
		 0.0000028249191933610997732134884831732080,
		 0.0000052619966475500025806047219888226911,
		 0.0000095660563720937891498760672094321504,
		 0.0000171390025138215001011925741410298063,
		 0.0000303721353159007172786038744071216398,
		 0.0000533144330199730628410098153685536815,
		 0.0000927727305168396906702150201057873596,
		 0.0001601089695752882413498063818822758630,
		 0.0002741571631011130146831419608340410510,
		 0.0004659284936908316930128467917102824680,
		 0.0007861555388635028741309694502552296580,
		 0.0013173193525957165981990005295537704555,
		 0.0021927170213724243173581210115258288695,
		 0.0036265296748535230252248606763743055126,
		 0.0059609930028473755913377739659719800613,
		 0.0097400298129714036937578519603775362157,
		 0.0158236368205604346909168636536868746134,
		 0.0255647567178983044495056784622954992869,
		 0.0410815423918531026012130723951143451700,
		 0.0656746545943673863224049509768054377901,
		 0.1044642348460616351570376912627491350349,
		 0.1653585520401501045929587227489854228679,
		 0.2605224824709100300400982375137459712278,
		 0.4086001179322432018067596387611573049981,
		 0.6380846956001181719268962322377802820483,
		 0.9924778169025332579173716285314554852448,
		 1.5384103118403334514111696629257153290382,
		 2.3793008307260859670196478932169270592567,
		 3.6817340408947272592000976576542825569049,
		 5.7409136583061526177239330959167773471563,
		 9.2321686012605892112867222998318084137281;
		
	for(int i = 0; i< d; i++){
		for(int j = 1; j<= 35; j++){
			DenseVector<Array<double> > tmp(maxvec[i]-minvec[i]+1);
			for(int l= tmp.firstIndex(); l<=tmp.lastIndex(); ++l){
				//std::cout << "i= " << i << "  j= " << j << "    " << (minvec(i)+l-tmp.firstIndex()) << endl;
				if(i == 1){
					tmp(l) = weights(j)*exp(-(minvec[i]+l-tmp.firstIndex())*(minvec[i]+l-tmp.firstIndex())*alpha(j));
				} else {
					tmp(l) = exp(-(minvec[i]+l-tmp.firstIndex())*(minvec[i]+l-tmp.firstIndex())*alpha(j));
				
				}

			}
			list.add(tmp);
		}
	}
	
};





double expsumApproximation(DimensionIndex list,DimensionIndex number_per_dim){
	DenseVector<Array<double> > weights(40);
	weights =  0.0004541105141216453968248119678948609979,
			   0.0078530680104171971521949987930788505963,
			 0.0358442951018547386594730302239453223478,
			 0.0800511636731031012202772917385473760987,
			 0.1168715918274595083080672451580905857327,
			 0.1317468847713531808620935140141661179314 ,
			 0.1268550732122552672878103285913908848670 ,
			 0.1109657824442714362255012552904709366430,
			 0.0915407675164530429816095620509930341768,
			 0.0728296028697006614068595181654686143702,
			 0.0566262616355938613997447837050680519155,
			 0.0433600894479567821162072083218674833915,
			 0.0328421513638694122585743456754769908912,
			 0.0246657172704435717145978130063199706967,
			 0.0183917192322462027550102406370235996746 ,
			 0.0136227224770988644418135728175744958435,
			 0.0100250532523774125004520276632780095838,
			 0.0073290867029904425973310905066682163778,
			 0.0053216028826095790999957683888726966970 ,
			 0.0038362256666861645005019393166179053622,
			 0.0027443224789343381565523340582385403152,
			 0.0019471510649087046089954599194064843726,
			 0.0013693829522368219844082718491921901727,
			 0.0009538877370448752111114987457026609396,
			 0.0006575902435281989381065818672453260341,
			 0.0004482091681493420014885365430801211772 ,
			 0.0003017066922888354188646730291625125409,
			 0.0002003052463432384262700059745642040820,
			 0.0001309533404038651058759093782370673542,
			 0.0000841450289679689093394985195683861284,
			 0.0000530163906346343672674381358229543301,
			 0.0000326579433267181342724703319572348084,
			 0.0000195946088007008467557055835143813537,
			 0.0000113949728253805224946867206361372001,
			 0.0000063797866387841876792628796961806114,
			 0.0000034061919433246302751298561507290710,
			 0.0000017093587076932593639466942818314337,
			 0.0000007873796897325934882231604846542004,
			 0.0000003186056947254193096300954077315363,
			 0.0000001024294180134919968361669149146366;
	DenseVector<Array<double> > alpha(40);
	alpha = 0.0328758165954420368502493336737568796480,
			 0.0570473982762439551973880462465427143570,
			 0.0966283152617544954693135257628444634292,
			 0.1636063258412157908437819728297846211262,
			 0.2786723770321616242787599732055880963344,
			 0.4784495396399847536917683543666868217770,
			 0.8285231470543625377543975329963643616793,
			 1.4475289223590702735247615073888027836801,
			 2.5522100303561091739872057626570267530042,
			 4.5426301651741418010516104786233881895896,
			 8.1651768488703318207497239278325196210062 ,
			14.8281159280173661825508180100996469263919,
			27.2202001694938288487574107321620431321207,
			50.5395968142454580715461354145645600510761,
			94.9698324334290623002785558526284148683771,
			180.7430429547675229534631213823558937292546,
			348.6595709978226785563659717581685981713235,
			682.3160325424825219386981700608885148540139 ,
			1355.9285760062651674529021761372860055416822,
			2739.2321520142396620300928589131217449903488,
			5632.4067005890412889179685862472979351878166,
			11804.1327408178242155400994306546635925769806,
			25254.1956632226066332691516436170786619186401,
			55255.5977602627446430005875299684703350067139,
			123898.7837399924441115217632614076137542724609,
			285403.3681506335250048778107156977057456970215 ,
			677314.9416636431666347561986185610294342041016,
			1661627.5164404432468927552690729498863220214844,
			4231169.6896806240947626065462827682495117187500,
			11239234.2676030291650022263638675212860107421875,
			31336866.8801911466143792495131492614746093750000,
			92434878.3605544002202805131673812866210937500000,
			291426049.6279058846412226557731628417968750000000,
			995657663.5402839354937896132469177246093750000000 ,
			3757776658.0723961133044213056564331054687500000000,
			16116802909.0660596201196312904357910156250000000000,
			82144510740.0004426687955856323242187500000000000000,
			538206356175.0550399124622344970703125000000000000000,
			5341579170601.2921476364135742187500000000000000000000,
			130450953442068.3405914306640625000000000000000000000000;

			double approx = 0.0;

			for(int i = 1; i<= 40; ++i){
				double prod = 0.0;
				for(int j = 0; j < list.length(); ++j){
					prod -= (list[j]*list[j]/(number_per_dim[j]*number_per_dim[j] + 0.0)) *alpha(i);
					//prod*= exp(-list(j)*list(j)*alpha(i)/list.length());
				}
				approx += weights(i)*exp(prod);
			}

			return approx;


};


double expsumApproximation_32(DimensionIndex list){
	DimensionIndex numpdim(list.length());
	for(int i = 0; i< numpdim.length(); ++i){
		numpdim[i] = 32;
	}
	return expsumApproximation(list,numpdim);
	
};

double expsumApproximation_64(DimensionIndex list){
	DimensionIndex numpdim(list.length());
	for(int i = 0; i< numpdim.length(); ++i){
		numpdim[i] = 64;
	}
	return expsumApproximation(list,numpdim);
	
};
double expsumApproximation_128(DimensionIndex list){
	DimensionIndex numpdim(list.length());
	for(int i = 0; i< numpdim.length(); ++i){
		numpdim[i] = 128;
	}
	return expsumApproximation(list,numpdim);
};
double expsumApproximation_256(DimensionIndex list){
	DimensionIndex numpdim(list.length());
	for(int i = 0; i< numpdim.length(); ++i){
		numpdim[i] = 256;
	}
	return expsumApproximation(list,numpdim);
};

double expsum(DimensionIndex list, DimensionIndex number_per_dim){
	double sum = 0.0;
	for(int i = 0; i< list.length(); ++i){
		sum += (list[i]/(number_per_dim[i] + 0.0))*(list[i]/(number_per_dim[i] + 0.0));
	}


	return exp(-sqrt(sum));
};


double expsum_32(DimensionIndex list){
	DimensionIndex maxval(list.length());
	for(int i = 0; i< maxval.length(); ++i){
		maxval[i] = 32;
	}
	return expsum(list,maxval);
};

double expsum_64(DimensionIndex list){
	DimensionIndex maxval(list.length());
	for(int i = 0; i< maxval.length(); ++i){
		maxval[i] = 64;
	}
	return expsum(list,maxval);
};

double expsum_128(DimensionIndex list){
	DimensionIndex maxval(list.length());
	for(int i = 0; i< maxval.length(); ++i){
		maxval[i] = 128;
	}
	return expsum(list,maxval);
};

double expsum_256(DimensionIndex list){
	DimensionIndex maxval(list.length());
	for(int i = 0; i< maxval.length(); ++i){
		maxval[i] = 256;
	}
	return expsum(list,maxval);
};

double inverseNormApproximation(DimensionIndex list){
DenseVector<Array<double> > weights(35);
	weights =  0.0003476617162031463656887349288667626857,
		 0.0003656677800731512040516686207159910427,
		 0.0004075639038080694722811763840587123697,
		 0.0004843102689961686493529664332240158853,
		 0.0006075484362407267528193426119921538575,
		 0.0007860621238313433964849066446726327806,
		 0.0010286920148215165835244833104043793437,
		 0.0013484039992954368087911095805929428959,
		 0.0017637932498565288986423389163381292288,
		 0.0022995948323400161628149234675083656398,
		 0.0029874656853097934646253271400578865880,
		 0.0038672614791946895807455515312045513099,
		 0.0049887848083617705673121177356588162866,
		 0.0064140049762051946081513006427121670328,
		 0.0082197972649914829865058925849330417890,
		 0.0105012852751240196090380674132092053696,
		 0.0133758979709688672253838190480401681981,
		 0.0169882797026389577249398925989920527968,
		 0.0215162197131539093567588515412691840822,
		 0.0271778000398019895103514419545942537582,
		 0.0342399988121945339294440423250076577233,
		 0.0430290321349169279780387271233665913428,
		 0.0539427787599617519049503139294543530013,
		 0.0674657211055092835482404310665227598065,
		 0.0841869982556718884943842105827993549383,
		 0.1048225290121168909272938649102080432840,
		 0.1302430953236043887537645713514677936473,
		 0.1615128713218071510428555659366622876405,
		 0.1999503386727134369738609689926001067306,
		 0.2472453213059679058203037954610792326093,
		 0.3057312880122441022998601434235155238639,
		 0.3791210128165055750809782353361043760742,
		 0.4747822440254026090267897020069653990504,
		 0.6124014951064692475124137671649293679366,
		 0.8770631816903148868821808725737554368607;
	DenseVector<Array<double> > alpha(35);
	alpha = 0.0000000235400856956305055878655722424688,
		 0.0000002190687866559873968066208739485593 ,
		 0.0000006537664626626933778312767948510514,
		 0.0000014417390459603761640437794711739619,
		 0.0000028249191933610997732134884831732080,
		 0.0000052619966475500025806047219888226911,
		 0.0000095660563720937891498760672094321504,
		 0.0000171390025138215001011925741410298063,
		 0.0000303721353159007172786038744071216398,
		 0.0000533144330199730628410098153685536815,
		 0.0000927727305168396906702150201057873596,
		 0.0001601089695752882413498063818822758630,
		 0.0002741571631011130146831419608340410510,
		 0.0004659284936908316930128467917102824680,
		 0.0007861555388635028741309694502552296580,
		 0.0013173193525957165981990005295537704555,
		 0.0021927170213724243173581210115258288695,
		 0.0036265296748535230252248606763743055126,
		 0.0059609930028473755913377739659719800613,
		 0.0097400298129714036937578519603775362157,
		 0.0158236368205604346909168636536868746134,
		 0.0255647567178983044495056784622954992869,
		 0.0410815423918531026012130723951143451700,
		 0.0656746545943673863224049509768054377901,
		 0.1044642348460616351570376912627491350349,
		 0.1653585520401501045929587227489854228679,
		 0.2605224824709100300400982375137459712278,
		 0.4086001179322432018067596387611573049981,
		 0.6380846956001181719268962322377802820483,
		 0.9924778169025332579173716285314554852448,
		 1.5384103118403334514111696629257153290382,
		 2.3793008307260859670196478932169270592567,
		 3.6817340408947272592000976576542825569049,
		 5.7409136583061526177239330959167773471563,
		 9.2321686012605892112867222998318084137281;

	double approx = 0.0;

	for(int i = 1; i<= 35; ++i){
		double prod = 0.0;
		for(int j = 0; j < list.length(); ++j){
			prod -= list[j]*list[j]*alpha(i);
			//prod*= exp(-list(j)*list(j)*alpha(i)/list.length());
		}
		approx += weights(i)*exp(prod);
	}

	return approx;
}; 


#endif // LAWA_HTUCKER_TENSORFUNCTIONS_H