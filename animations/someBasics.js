/*
 * Put <script src="egdod.js" data-scriptid="YOURINITSCRIPTID"></script> before your createCindy.
 */
(function(){
	var code = document.createTextNode(`

	computeDependentPoints() := (
		trafo = map([0,0], [3,0], [3,3], [0,3], a,b,c,d);
	
		vx = meet(join(a,b), join(c,d)).xy;
		vy = meet(join(a,d), join(b,c)).xy;
		m = meet(join(a,c), join(b,d)).xy;
	
		ad33 = (trafo * [0,1,1]).xy;
		ad50 = (trafo * [0,1.5,1]).xy;
		ad67 = (trafo * [0,2,1]).xy;
	
		bc33 = (trafo * [3,1,1]).xy;
		bc50 = (trafo * [3,1.5,1]).xy;
		bc67 = (trafo * [3,2,1]).xy;
		
		ab33 = (trafo * [1,0,1]).xy;
		ab50 = (trafo * [1.5,0,1]).xy;
		ab67 = (trafo * [2,0,1]).xy;
		
		dc33 = (trafo * [1,3,1]).xy;
		dc50 = (trafo * [1.5,3,1]).xy;
		dc67 = (trafo * [2,3,1]).xy;	
	);

	// *****************************************************************************

	a = [-22, -17];
	b = [14/3, -41/3];
	c = [2, -1];
	d = [-10, 4];

	computeDependentPoints();

	lineSizePrimary = 9;
	lineSizeSecondary = 4;
	
	pointSizePrimary = 0.5;
	pointSizeSecondary = 0.3;

	markerSizePrimary = 2.5 * pointSizePrimary;

	overshootPrimary = 0.1;
`);


	var scriptId = document.currentScript.dataset.scriptid;
	if (!scriptId) scriptId = 'csinit';

	var scriptElement = document.getElementById(scriptId);
	if (!scriptElement) {
		scriptElement = document.createElement("script");
		scriptElement.id = scriptId;
		scriptElement.type = "text/x-cindyscript";
		document.head.appendChild(scriptElement);
	}
	if (scriptElement.firstChild) {
		scriptElement.insertBefore(code, scriptElement.firstChild);
	} else {
		scriptElement.appendChild(code);
	}

})();
