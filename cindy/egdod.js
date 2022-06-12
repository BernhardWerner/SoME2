/*
 * Put <script src="egdod.js" data-scriptid="YOURINITSCRIPTID"></script> before your createCindy.
 */
(function(){
	var code = document.createTextNode(`
	canvasCorners = apply(screenbounds(), #.xy); //LO, RO, RU, LU
	canvasCorners = {
		"tl": canvasCorners_1,
		"tr": canvasCorners_2,
		"br": canvasCorners_3,
		"bl": canvasCorners_4
	};
	canvasCenter  = 0.5 * canvasCorners.tl + 0.5 * canvasCorners.br;
	canvasWidth   = dist(canvasCorners.tl, canvasCorners.tr);
	canvasHeight  = dist(canvasCorners.tl, canvasCorners.bl);
	[canvasLeft, canvasTop] = canvasCorners.tl;
	[canvasRight, canvasBottom] = canvasCorners.br;
	screenMouse() := [(mouse().x - canvasLeft) / canvasWidth, (mouse().y - canvasBottom) / canvasHeight];

	
	
	strokeSampleRateEBOW = 64;
	lineSampleRateEBOW = 64;
	texDelimitersEBOW = ["@", "€"];


	// ************************************************************************************************
	// Reverts a string.
	// ************************************************************************************************
	reverseString(string) := sum(tokenize(string, ""));


	// ************************************************************************************************
	// Gets the current time from the computer clock converted to seconds.
	// ************************************************************************************************
	computerSeconds() := (
		regional(actualTime);

		actualTime = time();

		actualTime_1 * 3600 + actualTime_2 * 60 + actualTime_3 + actualTime_4 * 0.001;
	);

	timeBufferEBOW = 0;
	scriptStartTimeEBOW = 0;
	// ************************************************************************************************
	// Sets up time-keeping variables. Will be automatically called when included.
	// Has to be called together with playanimation()!
	// ************************************************************************************************
	setupTime() := (
		timeBufferEBOW = computerSeconds();
		scriptStartTimeEBOW = timeBufferEBOW;
	);
	now() := computerSeconds() - scriptStartTimeEBOW;

	// ************************************************************************************************
	// Returns the duration ofthe last frame/tick in seconds.
	// Needs to run on every frame!
	// ************************************************************************************************
	deltaTime() := (
		regional(result);

		result = computerSeconds() - timeBufferEBOW;
		timeBufferEBOW = computerSeconds();

		result;
	);




	// ************************************************************************************************
	// Easing functions.
	// ************************************************************************************************
	easeInSine(x)        := 1 - cos((x * pi) / 2);
	easeOutSine(x)       := sin((x * pi) / 2);
	easeInOutSine(x)     := -(cos(pi * x) - 1) / 2;
   
	easeInQuad(x)        := x^2;
	easeOutQuad(x)       := 1 - (1 - x)^2;
	easeInOutQuad(x)     := if(x < 0.5, 2 * x^2, 1 - (-2 * x + 2)^2 / 2);
	   
	easeInCubic(x)       := x^3;
	easeOutCubic(x)      := 1 - (1 - x)^3;
	easeInOutCubic(x)    := if( x < 0.5, 4 * x^3, 1 - (-2 * x + 2)^3 / 2);
	   
	easeInQuart(x)       := x^4;
	easeOutQuart(x)      := 1 - (1 - x)^4;
	easeInOutQuart(x)    := if(x < 0.5, 8 * x^4, 1 - (-2 * x + 2)^4 / 2);
	   
	easeInQuint(x)       := x^5;
	easeOutQuint(x)      := 1 - (1 - x)^5;
	easeInOutQuint(x)    := if(x < 0.5, 16 * x^5, 1 - (-2 * x + 2)^5 / 2);
   
	easeInExpo(x)        := if(x == 0, 0, 2^(10 * x - 10));
	easeOutExpo(x)       := if(x == 1, 1, 1 - 2^(-10 * x));
	easeInOutExpo(x)     := if(x == 0, 0, if(x == 1, 1, if(x < 0.5, 2^(20 * x - 10) / 2, (2 - 2^(-20 * x + 10)) / 2)));
   
	easeInCirc(x)        := 1 - sqrt(1 - x^2);
	easeOutCirc(x)       := sqrt(1 - (x - 1)^2);
	easeInOutCirc(x)     := if(x < 0.5, (1 - sqrt(1 - 4 * x^2)) / 2, (sqrt(1 - (-2 * x + 2)^2) + 1) / 2);

	easeInBack(x)        := 2.70158 * x^3 - 1.70158 * x^2;
	easeOutBack(x)       := 1 - easeInBack(1 - x);
	easeInOutBack(x)     := if(x < 0.5, 4 * x^2 * ((1.70158 * 1.525 + 1) * 2 * x - 1.70158 * 1.525) / 2, ((2 * x - 2)^2 * ((1.70158 * 1.525 + 1) * (2 * x - 2) + 1.70158 * 1.525) + 2) / 2);

	easeInElastic(x)     := if(x == 0, 0, if(x == 1, 1, -2^(10 * x - 10) * sin(2 * pi / 3 * (10 * x - 10.75))));
	easeOutElastic(x)    := 1 - easeInElastic(1 - x);
	easeInOutElastic(x)  := if(x == 0, 0, if(x == 1, 1, if(x < 0.5, -2^(20 * x - 10) * sin(4 * pi / 9 * (20 * x - 11.125)) / 2, 2^(-20 * x + 10) * sin(4 * pi / 9 * (20 * x - 11.125)) / 2 + 1)));



	// ************************************************************************************************
	// Basic animation functionlity.
	// ************************************************************************************************

	setupAnimationTrack(s, e) := {
		"start":    s,
		"end":      e,
		"duration": e - s,
		"timeLeft": e - s,
		"progress": 0,
		"running":  true,
		"looping":  false
	}; 

	trackStarted(track, delay) := now() >= track.start + delay;
	trackEnded(track, delay) := now() > track.end + delay;
	trackRunning(track, delay) := and(trackStarted(track, delay), not(trackEnded(track, delay)));
	trackStarted(track) := trackStarted(track, 0);	
	trackEnded(track) := trackEnded(track, 0);
	trackRunning(track) := trackRunning(track, 0);


	// Needs to run on every frame!
	updateAnimationTrack(track, delta) := (
		if(track.running & trackStarted(track),
			track.timeLeft = track.timeLeft - delta;	
			track.progress = 1 - track.timeLeft / track.duration;
			if(track.timeLeft <= 0,
				if(track.looping,
					track.timeLeft = track.end - track.start;
				, // else //
					track.timeLeft = 0;
					track.progress = 1;
					track.running = false;		
				);
			);
		);
	);

	tween(obj, prop, from, to, track, easing) := (
		regional(t);

		t = track.progress;
		
		if(trackStarted(track),
			if(t < 1,
				if(easing != "none",
					t = parse(easing + "(" + t + ")");
				);
			
				if(contains(keys(obj), prop),
					obj_prop = lerp(from, to, t);
				);	
			,if(t >= 1,
				if(contains(keys(obj), prop),
					obj_prop = to;
				);	
			));
		);
	);
	tween(obj, prop, from, to, track) := tween(obj, prop, from, to, track, "none");

	
	tweenMany(list, prop, from, to, track, delay, easing) := (
		regional(t);


			forall(1..length(list),
				t = 1 - (track.timeLeft + (# - 1) * delay) / track.duration;
				
				if(trackStarted(track),
					if(track.timeLeft < track.duration - (# - 1) * delay,
						if(track.timeLeft > - (# - 1) * delay,
							if(easing != "none",
								t = parse(easing + "(" + t + ")");
							);
						
							if(contains(keys(list_#), prop),
								list_#_prop = lerp(from, to, t);
							);	
						,if(track.timeLeft <= - (# - 1) * delay,
							if(contains(keys(list_#), prop),
								list_#_prop = to;
							);	
						));
					);
				);
			);

	);


	tweenMany(list, prop, from, to, track, delay) := tweenMany(list, prop, from, to, track, delay, "none");

	// ************************************************************************************************
	// Setting up several animation tracks with delay.
	// ************************************************************************************************
	setupAnimationCascade(number, start, duration, step) := (
		apply(0..number - 1,
			setupAnimationTrack(start + # * step, start + # * step + duration);	
		);
	);

	// ************************************************************************************************
	// Updating all tracks in a cascade. Has to run on every frame.
	// ************************************************************************************************
	updateAnimationCascade(cascade, delta) := (
		forall(cascade, track,
			updateAnimationTrack(track, delta);
		);
	);

	// ************************************************************************************************
	// Updating all tracks in a cascade. Has to run on every frame.
	// ************************************************************************************************
	runAnimationCascade(cascade) := forall(cascade, track, track.running = true);




	// ************************************************************************************************
	// Rendering various animation objects.
	// ************************************************************************************************
	// Text objects needs
	// pos
	// offset
	// text
	// percentVisible
	// size
	// alpha
	// color 
	// fontFamily
	// align
	// ************************************************************************************************
	drawTextObject(obj) := (
		drawtext(obj.pos + obj.offset, substring(obj.text, 0, round(obj.percentVisible * length(obj.text))), size->obj.size, color->obj.color, align->obj.align, family->obj.fontFamily, align->obj.align, alpha->obj.alpha);
	);


	createTextObject(pos, text, size) := {
		"pos": pos,
		"offset": [0,0],
		"text": text,
		"percentVisible": 1,
		"size": size,
		"alpha": 1,
		"color ": (0,0,0),
		"align": "left"
	};





















	slideIn(obj, dir, dist, track) := (
		if(dir != [0,0],
			dir = dir / abs(dir);
			tween(obj, "offset", dist * dir, [0,0], track, "easeInOutBack");	
		);
	);
	slideOut(obj, dir, dist, track) := (
		if(dir != [0,0],
			dir = dir / abs(dir);
			tween(obj, "offset", [0,0], dist * dir, track, "easeInOutBack");	
		);
	);

	fadeIn(obj, dir, dist, track) := (
		slideIn(obj, dir, dist, track);
		tween(obj, "alpha", 0, 1, track, "easeInQuint");
	);
	fadeOut(obj, dir, dist, track) := (
		slideOut(obj, dir, dist, track);
		tween(obj, "alpha", 1, 0, track, "easeOutQuint");
	);






	// ************************************************************************************************
	// Stroke object needs
	// pos
	// stroke (list of points, relative to pos)
	// drawStart
	// drawEnd
	// scale
	// lineColor
	// lineAlpha
	// lineSize
	// fillColor
	// fillAlpha
	// arrow
	// arrowSize
	// ************************************************************************************************
	drawStrokeObject(obj) := (
		regional(absoluteStroke, ratio, n);
		
		if(obj.scale ~!= 0,
			absoluteStroke = apply(obj.stroke, obj.pos + obj.offset + obj.scale * #);
			n = length(absoluteStroke);

			ratio = max(1,floor(obj.drawStart * n))..min(n, ceil(obj.drawEnd * n));
			fillpoly(absoluteStroke_ratio, color->obj.fillColor, alpha->obj.fillAlpha);
			connect(absoluteStroke_ratio, size->obj.scale * obj.lineSize, color->obj.lineColor, alpha->obj.lineAlpha);
			if(obj.arrow, 
				if(length(ratio) >= 3,
					connect(arrowTip(absoluteStroke_(ratio_(-1)), absoluteStroke_(ratio_(-1)) - absoluteStroke_(ratio_(-3)), obj.scale * obj.arrowSize), size->obj.scale * obj.lineSize, color->obj.lineColor, alpha->obj.lineAlpha);
				,if(length(ratio) >= 2,
					connect(arrowTip(absoluteStroke_(ratio_(-1)), absoluteStroke_(ratio_(-1)) - absoluteStroke_(ratio_(-2)), obj.scale * obj.arrowSize), size->obj.scale * obj.lineSize, color->obj.lineColor, alpha->obj.lineAlpha);
				));
			);
		);
	);
	arrowTipAngleEBOW = pi/ 6;
	arrowTip(tipPos, dir, size) := (
		if(abs(dir) > 0, dir = dir / abs(dir));

		[
			tipPos - size * rotation(arrowTipAngleEBOW) * dir,
			tipPos,
			tipPos - size * rotation(-arrowTipAngleEBOW) * dir
		];		
	);


	// ************************************************************************************************
	// Flipbook object needs
	// pos
	// offset
	// scale
	// flipbook
	// index
	// flipx
	// flipy
	// rotation
	// alpha
	// ************************************************************************************************
	drawFlipbookObject(obj) := (
		drawimage(obj.pos + obj.offset, obj.flipbook_(obj.index), scale->obj.scale, rotation->obj.rotation, flipx->obj.flipx, flipy->obj.flipy, alpha->obj.alpha);
	);
	

	// ************************************************************************************************
	// Cycles through flipbook of a flipbook object. Flipbook pages are equdistantly 
	// spread across animation track as key frames.
	// ************************************************************************************************
	animateFlipbook(obj, track) := (
		regional(n);
		n = length(obj.flipbook);

		obj.index = floor(lerp(1, n + 0.9999, track.timeLeft, track.end, track.start));
	);


	
	
	// ************************************************************************************************
	// Setting up a stroke object.
	// ************************************************************************************************
	createRootStrokeObject(pos, lineColor, fillColor, fillAlpha) := {
		"pos": pos,
		"offset": [0,0],
		"stroke": apply(1..strokeSampleRateEBOW, [0,0]),
		"drawStart": 0,
		"drawEnd": 0,
		"scale": 1,
		"lineSize": 0,
		"lineColor": lineColor,
		"lineAlpha": 1,
		"fillColor": fillColor,
		"fillAlpha": fillAlpha,
		"arrow": false,
		"arrowSize": 0
	};
	createRootStrokeObject(pos, lineColor) := createRootStrokeObject(pos, lineColor, (1,1,1), 0);
	
	// ************************************************************************************************
	// Zero strokes.
	// ************************************************************************************************
	zeroStrokeCenter() := apply(1..strokeSampleRateEBOW, [0,0]);
	zeroStrokeRight() := apply(1..strokeSampleRateEBOW, [1,0]);
	
	// ************************************************************************************************
	// Creates stroke around a circle.
	// ************************************************************************************************
	sampleCircle(rad, angle) := apply(0..strokeSampleRateEBOW - 1, rad * [cos(angle * # / (strokeSampleRateEBOW - 1)), sin(angle * # / (strokeSampleRateEBOW - 1))]);
	
	
	// ************************************************************************************************
	// Subdivides the distance between two points.
	// ************************************************************************************************
	subdivideSegment(p, q, n) := apply(1..n, lerp(p, q, #, 1, n));
	
	// ************************************************************************************************
	// Creates stroke around a polygon.
	// ************************************************************************************************
	samplePolygonFREE(poly, nop, closed) := (
		regional(pairs, dists, totalDist, effectiveNumber, splitNumbers, stepSize, index, sr);
		
		if(closed, 
			poly = poly :> poly_1;
		);
		
		sr = if(length(poly) == 2, lineSampleRateEBOW, nop);
		
		pairs = consecutive(poly);
		
		dists = apply(pairs, dist(#_1, #_2));
		totalDist = sum(dists);
		
		effectiveNumber = sr - length(poly);

		splitNumbers = apply(dists, floor((sr - 1) * # / totalDist));


		
		if(sum(splitNumbers) < effectiveNumber,
			forall(1..effectiveNumber - sum(splitNumbers),
				index = randchoose(1..length(splitNumbers));
				splitNumbers_index = splitNumbers_index + 1;
			);
		);
		if(sum(splitNumbers) > effectiveNumber,
			forall(1..sum(splitNumbers) - effectiveNumber,
				index = randchoose(1..length(splitNumbers));
				splitNumbers_index = splitNumbers_index - 1;
			);
		);


		flatten(apply(1..length(pairs), pop(subDivideSegment(pairs_#_1, pairs_#_2, splitNumbers_# + 2)) )) :> poly_(-1);
	
	);
	samplePolygonFREE(poly, nop) :=	samplePolygonFREE(poly, nop, true);

	samplePolygon(poly) := samplePolygonFREE(poly, strokeSampleRateEBOW);
	samplePolygon(poly, closed) := samplePolygonFREE(poly, strokeSampleRateEBOW, closed);

	




	// ************************************************************************************************
	// Resampling via centripetal Catmull-Rom splines.
	// ************************************************************************************************
	resample(stroke, nop) := sampleCatmullRomSplineFREE(stroke, nop);
	resample(stroke) := sampleCatmullRomSplineFREE(stroke, strokeSampleRateEBOW);




	// ************************************************************************************************
	// Creates stroke as Bezier curves.
	// ************************************************************************************************
	bezier(controls, t) := (
		regional(n);

		n = length(controls);

		sum(apply(1..n, binom(n - 1, # - 1) * t^(# - 1) * (1 - t)^(n - #) * controls_#));
	);

	sampleBezierCurve(controls) := (
		regional(t, sr);

		sr = if(length(controls) == 2, lineSampleRateEBOW, strokeSampleRateEBOW);
		apply(0..sr - 1, 
			t = # / (sr - 1);

			bezier(controls, t);
		);
	);

	// ************************************************************************************************
	// Discrete derivative of a discrete curve
	// ************************************************************************************************
	derive(curve) := apply(consecutive(curve), #_2 - #_1);


	// ************************************************************************************************
	// Creates stroke as Catmull-Rom curves.
	// ************************************************************************************************
	catmullRom(controls, alpha, t) := (
		regional(a, b, c, p, q, knots);

		knots = [0];
		forall(2..4,
			knots = knots :> knots_(-1)	+ dist(controls_#, controls_(# - 1))^alpha;
		);
		t = lerp(knots_2, knots_3, t);

		a = lerp(controls_1, controls_2, t, knots_1, knots_2);
		b = lerp(controls_2, controls_3, t, knots_2, knots_3);
		c = lerp(controls_3, controls_4, t, knots_3, knots_4);

		p = lerp(a, b, t, knots_1, knots_3);
		q = lerp(b, c, t, knots_2, knots_4);

		lerp(p, q, t, knots_2, knots_3);
	);


	
	
	sampleCatmullRomCurve(controls, alpha) := (
		regional(t);
		apply(0..strokeSampleRateEBOW - 1, 
			t = # / (strokeSampleRateEBOW - 1);

			catmullRom(controls, alpha, t);
		);
	);
	sampleCatmullRomCurve(controls) := sampleCatmullRomCurve(controls, 0.5);
			
	sampleCatmullRomSplineGeneralFREE(points, alpha, nop) := (
		regional(dists, traj, before, after, cutTimes, piece, controls, t);

		dists    = apply(derive(points), abs(#));
		traj     = sum(dists);
		before   = 2 * points_1    - points_2;
		after    = 2 * points_(-1) - points_(-2);
		cutTimes = 0 <: apply(1..(length(dists) - 1), sum(dists_(1..#))) / traj;
	  
		apply(0..(nop - 1), i,
		  piece = select(1..(length(points) - 1), cutTimes_# * (nop - 1) <= i)_(-1);
	  
		
		  if(piece == 1,
			controls = [before, points_1, points_2, points_3];
			t = i / (nop - 1) * traj / dists_1;
		  ,if(piece == length(points) - 1,
			controls = [points_(-3), points_(-2), points_(-1), after];
			t = (i / (nop - 1) - cutTimes_(-1)) * traj / dists_(-1);
		  , // else //
			controls = [points_(piece - 1), points_(piece), points_(piece + 1), points_(piece + 2)];
			t = (i / (nop - 1) - cutTimes_piece) * traj / dists_piece;
		  ));

		  catmullRom(controls, alpha, t);
		);
	);
	sampleCatmullRomSplineGeneral(points, alpha) := sampleCatmullRomSplineGeneralFREE(points, alpha, strokeSampleRateEBOW);
	sampleCatmullRomSplineFREE(points, nop)      := sampleCatmullRomSplineGeneralFREE(points, 0.5, nop);
	sampleCatmullRomSpline(points)               := sampleCatmullRomSplineGeneralFREE(points, 0.5, strokeSampleRateEBOW);
	


	

	// ************************************************************************************************
	// Creates a stroke based on a function graph.
	// ************************************************************************************************
	sampleFunctionGraph(func, start, end) := (
		apply(0..strokeSampleRateEBOW-1,
			t = lerp(start, end, #, 0, strokeSampleRateEBOW-1);

			(t, parse(func + "(" + t + ")"));	
		);
	);
	sampleCurve(curve, start, end) := (
		apply(0..strokeSampleRateEBOW-1,
			t = lerp(start, end, #, 0, strokeSampleRateEBOW-1);

			parse(curve + "(" + t + ")");	
		);
	);


	// ************************************************************************************************
	// Draws a stroke object as a stroke.
	// ************************************************************************************************
	constructStrokeDraw(obj, lineSize, arrowSize, track) := (
		tween(obj, "lineSize", 0, lineSize, track, "easeOutCirc");
		tween(obj, "drawEnd", 0, 1, track, "easeInOutCubic");
		tween(obj, "arrowSize", 0, arrowSize, track);	
	);
	destroyStrokeDrawBackwards(obj, lineSize, arrowSize, track) := (
		tween(obj, "lineSize", lineSize, 0, track, "easeInCirc");
		tween(obj, "drawEnd", 1, 0, track, "easeInOutCubic");
		tween(obj, "arrowSize", arrowSize, 0, track);	
	);
	destroyStrokeDrawForwards(obj, lineSize, arrowSize, track) := (
		tween(obj, "lineSize", lineSize, 0, track, "easeInCirc");
		tween(obj, "drawStart", 0, 1, track, "easeInOutCubic");
		tween(obj, "arrowSize", arrowSize, 0, track);	
	);



	wormStrokeDraw(obj, lineSize, arrowSize, wormSize, track) := (
		if(track.progress <= 0.5,
			tween(obj, "lineSize", 0, lineSize, track, "easeOutCirc");
			tween(obj, "arrowSize", 0, arrowSize, track, "easeOutCirc");
		, // else //
			tween(obj, "lineSize", lineSize, 0, track, "easeInCirc");
			tween(obj, "arrowSize", arrowSize, 0, track, "easeInCirc");
		);
		tween(obj, "drawEnd", 0, 1 + wormSize, track, "easeInOutQuad");
		obj.drawStart = obj.drawEnd - wormSize;
			
	);

	
	constructStrokeDrawMany(list, lineSize, arrowSize, track, delay) := (
		tweenMany(list, "lineSize", 0, lineSize, track, delay, "easeOutCirc");
		tweenMany(list, "drawEnd", 0, 1, track, delay, "easeInOutCubic");
		tweenMany(list, "arrowSize", 0, arrowSize, track, delay);	
	);
	destroyStrokeDrawBackwardsMany(obj, lineSize, arrowSize, track, delay) := (
		tweenMany(obj, "lineSize", lineSize, 0, track, delay, "easeInCirc");
		tweenMany(obj, "drawEnd", 1, 0, track, delay, "easeInOutCubic");
		tweenMany(obj, "arrowSize", arrowSize, 0, track, delay);	
	);
	destroyStrokeDrawForwardsMany(obj, lineSize, arrowSize, track, delay) := (
		tweenMany(obj, "lineSize", lineSize, 0, track, delay, "easeInCirc");
		tweenMany(obj, "drawStart", 0, 1, track, delay, "easeInOutCubic");
		tweenMany(obj, "arrowSize", arrowSize, 0, track, delay);	
	);


	// ************************************************************************************************
	// Grows a stroke.
	// ************************************************************************************************
	constructStrokeGrow(obj, track) := (
		tween(obj, "scale", 0, 1, track, "easeOutBack");
	);

	constructStrokeGrowMany(list, track, delay) := (
		tweenMany(list, "scale", 0, 1, track, delay, "easeOutBack");
	);

	destroyStrokeShrink(obj, track) := (
		tween(obj, "scale", 1, 0, track, "easeInBack");
	);

	destroyStrokeShrinkMany(list, track, delay) := (
		tweenMany(list, "scale", 1, 0, track, delay, "easeInBack");
	);

		
	
	// ************************************************************************************************
	// Transformation matrices.
	// ************************************************************************************************
	rotation(alpha) := [[cos(alpha), -sin(alpha)], [sin(alpha), cos(alpha)]];

	
	rotate(point, alpha, center) := rotation(alpha) * (point - center) + center;
	rotate(point, alpha) := rotate(point, alpha, [0,0]);

	












    	
		// ************************************************************************************************
		// Takes two arrays of the same length and pairs elements with the same index together.
		// ************************************************************************************************
		zip(a, b) := transpose([a, b]);



		// *************************************************************************************************
		// Removes the first i elements of a list.
		// *************************************************************************************************
		bite(list, i) := list_((i + 1)..length(list));
		bite(list) := bite(list, 1);


		// *************************************************************************************************
		// Intersecting and joining a list of lists.
		// *************************************************************************************************
		cap(list) := (
		    regional(res);
		    if(list == [],
		        [];
		    , // else //
		        res = list_1;
		        forall(bite(list),
		            res = res ~~ #;
		        );

		        res;
		    );
		);
		cup(list) := set(flatten(list));



		// *************************************************************************************************
		// Generates all triples of consecutive elements of a list.
		// *************************************************************************************************
		consectriples(list) := (
		    regional(res);

		    res = [];
		    if(length(list) <= 2,
		        res = [];
		    , // else //
		        forall(1..(length(list) - 2),
		            res = res :> list_[#, # + 1, # + 2];
		        );
		    );
		);

		// *************************************************************************************************
		// Returns a list of length n where every entry is the object x.
		// *************************************************************************************************
		const(n, x) := if(n == 0, [], apply(1..n, x));


		
		// *************************************************************************************************
		// Finds first index at which x appears in list. Returns 0, when x is not in list.
		// *************************************************************************************************
		findin(list, x) := (
		    regional(occs);

		    occs = select(1..length(list), list_# == x);
		    if(length(occs) == 0, 0, occs_1);
		);

		// ************************************************************************************************
		// Returns the number of elements in list that are equal to x.
		// ************************************************************************************************
		frequency(list, x) := length(select(list, # == x));



		// *************************************************************************************************
		// Checks whether a list is constant or constant with a specific value.
		// *************************************************************************************************
		isconst(list) := (
		         list == const(length(list), list_1);
		);
		isconst(list, x) := (
		         list == const(length(list), x);
		);

		// *************************************************************************************************
		// For two lists consisting of distinct elements this function returns the permutation mapping the
		// first list to the second. It does it such that list1_result = list2.
		// Please make sure yourself, that the lists are compatible, i.e. that such a permutation exists.
		// *************************************************************************************************
		findperm(list1, list2) := (
		    apply(list2, e, select(1..length(list1), list1_# == e)_1);
		);

		// *************************************************************************************************
		// Removes last i elements of an array.
		// *************************************************************************************************
		pop(list) := list_(1..(length(list) - 1));
		pop(list, i) := list_(1..(length(list) - i));

		// *************************************************************************************************
		// Provides a list of all integers 1...n randomly sorted.
		// *************************************************************************************************
		randomindex(n) := randsort(1..n);

		// *************************************************************************************************
		// Sorts a list randomly.
		// *************************************************************************************************
		randsort(list) := (
			regional(l, temp, i);

			l = length(list);

			while(l > 0,
				i = randomint(l) + 1;

				temp = list_l;
				list_l = list_i;
				list_i = temp;
				l = l - 1;
			);

			list;
		);

		

		// *************************************************************************************************
		// Chooses randomly k elements of a list.
		// *************************************************************************************************
		randchoose(list, k) := (
				regional(res, i);

				if(k > length(list),
					randsort(list);
				, // else //
					res = [];
					forall(1..k,
						i = randomint(length(list)) + 1;
						res = res :> list_i;
						list = remove(list, i);
					);

					res;
				);
		);

		randchoose(list) := randchoose(list, 1)_1;



		// *************************************************************************************************
		// Removes the i-th element of an array.
		// *************************************************************************************************
		remove(arr, i) := (
		   if(i <= 1,
		    bite(arr);
		  ,if(i >= length(arr),
		    pop(arr);
		  , // else //
		    arr_((1..(i - 1)) ++ ((i + 1)..length(arr)));
		  ));
		);


		// *************************************************************************************************
		// Removes elements from a that lie in b.
		// *************************************************************************************************
		setMinus(a, b) := a -- b; //select(a, !contains(b, #));














		// ************************************************************************************************
		// Gets color value out of a gradient. Gradient time values have to be between 0 and 1.
		// Gradient has to be an array of JSONs with keys color and t.
		// ************************************************************************************************
		evalGradient(grad, time) := (
			regional(index);
			
			grad = sort(grad, #.t);

			time = clamp(time, 0, 1);

			index = 1;
			forall(1..(length(grad) - 1),
				if(time > grad_#.t, index = #);
			);

			lerp(grad_index.color, grad_(index + 1).color, time, grad_index.t, grad_(index + 1).t);
		);


		// ************************************************************************************************
		// Color transformation
		// ************************************************************************************************
		rgb2hsv(vec) := (
			regional(cMax, cMin, delta, maxIndex, h, s);

			maxIndex = 1;
			cMax = vec_1;
			forall(2..3, 
				if(vec_# > cMax,
					maxIndex = #;
					cMax = vec_#;
				);
			);

			cMin = min(vec);
			delta = cMax - cMin;

			h =  if(delta ~= 0,
					0;
				,if(maxIndex == 1,
					mod((vec_2 - vec_3) / delta, 6);
				,if(maxIndex == 2,
					2 + (vec_3 - vec_1) / delta
				,if(maxIndex == 3,
					4 + (vec_1 - vec_2) / delta
					
				)))) / 6;

			s = if(cMax ~= 0, 0, delta / cMax);

			[h, s, cMax];

		);

		hsv2rgb(vec) := (
			regional(c, x, m, r, g, b);

			vec_1 = vec_1 * 6;
			c = vec_2 * vec_3;
			x = c * (1 - abs(mod(vec_1, 2) - 1));
			m = vec_3 - c;

			[r, g, b] =  if(vec_1 < 1, 
							[c, x, 0];
						,if(vec_1 < 2, 
							[x, c, 0];
						,if(vec_1 < 3, 
							[0, c, x];
						,if(vec_1 < 4, 
							[0, x, c];
						,if(vec_1 < 5, 
							[x, 0, c];
						,if(vec_1 < 6, 
							[c, 0, x];
						))))));

			[r + m, g + m, b + m];
		);
		
		deca2hexa(digit) := ["0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F"]_(digit + 1);
		hexa2deca(digit) := (
			regional(x, y);

			x = findin(["0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F"], digit) - 1;
			y = findin(["0","1","2","3","4","5","6","7","8","9","a","b","c","d","e","f"], digit) - 1;

			if(x == -1, y, x);
		);

		rgb2hex(vec) := (
			regional(a, b);
			vec = round(255 * vec);

			sum(apply(vec,
				a = mod(#, 16);
				b = (# - a) / 16;
				deca2hexa(b) + deca2hexa(a);
			));			
		);

	

		hex2rgb(string) := (
			regional(digits);

			digits = tokenize(string, "");
			apply([1,3,5],
				16 * hexa2deca(text(digits_#)) + hexa2deca(text(digits_(# + 1)));
			) / 255;
		);






        // ************************************************************************************************
		// Checks whether float a lies between floats c and d.
		// ************************************************************************************************
		between(a, c, d) := (a >= c) & (a <= d);


		// ************************************************************************************************
		// Clamps x between the values a and b.
		// ************************************************************************************************
		clamp(x, a, b) := if(a <= b, min(max(x, a), b), min(max(x, b), a));

		// ************************************************************************************************
		// t-conorm dual to normal multiplication. (I.e. Hamacher p with p = 1).
		// ************************************************************************************************
		oplus(x, y) := x + y - x * y;


		pNorm(p, v) := (
			v = apply(v, abs(#));

			if(p == 0,
				max(v);
			, // else //
				sum(apply(v, #^p))^(1 / p);
			);
		);

		distPointSet(point, set) := (
			min(apply(set, dist(point,#)));
		);

		// *************************************************************************************************
		// Computes the binomial n over k.
		// *************************************************************************************************
		binom(n, k) := (
		    if((n < 0) % (k < 0) % (k > n),
		        err("binom: wrong numbers")
		    , // else //
		        faculty(n) / faculty(k) / faculty(n - k)
		    );
		);
		


		// *************************************************************************************************
		// Computes the convexhulls of a list of points
		// *************************************************************************************************
		cross(o, a, b) := (
			(a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)
		  );
		convexHull(points) := (
			regional(ordered, upper, lower);
		  
			ordered = set(sort(points));
			if(length(ordered) <= 3,
			  ordered;
			, // else //
			  lower = [];
			  forall(ordered,
				while((length(lower) > 1) & (cross(lower_(-2), lower_(-1), #) <= 0),
				  lower = pop(lower);
				);
				lower = lower :> #;
			  );
			  upper = [];
			  forall(reverse(ordered),
				while((length(upper) > 1) & (cross(upper_(-2), upper_(-1), #) <= 0),
				  upper = pop(upper);
				);
				upper = upper :> #;
			  );
		  
			  pop(lower) ++ pop(upper);
			);
		  );
		



		// *************************************************************************************************
		// Area of a (simple?) polygon.
		// *************************************************************************************************
		areaOfPolygon(points) := (
			0.5 * abs(sum(apply(1..(length(points) - 1),
			  (points_#).x * (points_(# + 1)).y - (points_(# + 1)).x * (points_#).y;
			)));
		  );


		// *************************************************************************************************
		// Lagrange interpolation for a list of points.
		// *************************************************************************************************
		lagrange(list, x) := sum(apply(list, p,
			p.y * product(apply(list -- [p], q, (x - q.x) / (p.x - q.x)));
		));

		
		sampleLagrangeInterpolationFREE(points, nop) := (
			regional(dists, traj, cutTimes, piece, t, start, end);

			[start, end] = [min(points).x, max(points).x];

			apply(1..nop, [lerp(start, end, #, 1, nop), lagrange(points, lerp(start, end, #, 1, nop))] );
		);
		sampleLagrangeInterpolation(points) := sampleLagrangeInterpolationFREE(points, strokeSampleRateEBOW);




		// *************************************************************************************************
		// Computes coefficients of polynomial interpolation via Newton basis / divided differences.
		// Output are coefficients for Newton basis ordered from lowest to highest degree.
		// *************************************************************************************************
		dividedDifferences(points) := (
			regional(recurMatrix, n);

			n = length(points);

			recurMatrix = zeromatrix(n, n);

			forall(1..n, i, forall(1..i, j,
				recurMatrix_i_j = if(j == 1,
					points_i.y;
				, // else //
					(recurMatrix_i_(j - 1) - recurMatrix_(i - 1)_(j - 1)) / (points_i.x - points_(i - j + 1).x);
				);
			));

			apply(1..n, recurMatrix_#_#);
		);

		// *************************************************************************************************
		// Evaluates polynomial given by list of coefficients of Newton basis at value x.
		// Coefficients have to be sorted from lowest to highest degree!
		// *************************************************************************************************
		newtonInterpolation(points, coeffs, x) := if(length(coeffs) == 1, coeffs_1, 
			newtonInterpolation(bite(points), bite(coeffs), x) * (x - points_1.x) + coeffs_1;
		);


		// *************************************************************************************************
		// Evaluates polynomial given by list of coefficients via Horner scheme at value x.
		// Coefficients have to be sorted from lowest to highest degree!
		// *************************************************************************************************
		horner(coeffs, x) := if(length(coeffs) == 1, coeffs_1, horner(bite(coeffs), x) * x + coeffs_1);


		// *************************************************************************************************
		// W.I.P. - DOESN'T WORK!!
		// The GJK collison detection algorithm for finite convex polygons in 2D.
		// *************************************************************************************************
		/*
		GJKcollision(shapeA, shapeB) := (
			regional(dir, simplex, found, a, b, c, perpB, perpC);

			if(shapeA == [] % shapeB == [],
				false;
			, // else //
				dir = sum(shapeB) / length(shapeB) - sum(shapeA) / length(shapeA);
				simplex = [polySupport(shapeA, shapeB, dir)];
				dir = (0,0) - simplex_1;

				result = false;
				searching = true;
				while(searching,
					a = polySupport(shapeA, shapeB, dir);
					if(a * dir < 0,
						searching = false;
					, // else //
						simplex = simplex :> a;
						if(length(simplex) == 2,
							[b, a] = simplex;
							if(triangleheight(a,b,(0,0)) ~= 0,
								result = true;
								searching = false;
							, // else // 
								dir = perp(b-a) * sign(triangleheight(a,b,(0,0)));
							);
						, // else //
							[c, b, a] = simplex;
							if(triangleheight(a,b,(0,0)) * triangleheight(a,c,(0,0)) * triangleheight(b,c,(0,0)) ~= 0,
								result = true;
								searching = false;
							, // else // 
								perpB = perp(b-a) * sign(triangleheight(a,b,(0,0)));
								perpC = perp(c-a) * sign(triangleheight(a,c,(0,0)));

								 if(perpB * a < 0,
									simplex = simplex_[2,3];
									dir = perpB;
								,if(perpC * a < 0,
									simplex = simplex_[1,3];
									dir = perpC;
								, // else //
									result = true;
									searching = false;
								));
							);
						);
					);
				);

				result;
			);
		);
		polySupport(shape, dir) := sort(shape, # * dir)_(-1);
		polySupport(shapeA, shapeB, dir) := polySupport(shapeA, dir) - polySupport(shapeB, -dir);
		*/

		lazyCollision(shapeA, shapeB) := (
			regional(mink);
			
			mink = convexHull(flatten(apply(shapeA, a, apply(shapeB, b, b - a))));

			result = true;
			forall(cycle(mink),
				result = result & (triangleheight(#_1, #_2, (0,0)) >= 0);
			);

			result;
		);

		// CONVEX POLYGONS ONLY!!!
		pointInPolygon(point, poly) := (
			regional(resultForwards, resultBackwards);

			resultForwards = true;
			resultBackwards = true;
			forall(cycle(poly),
				resultForwards  = and(resultForwards , triangleheight(#_1, #_2, point) >= 0);
				resultBackwards = and(resultBackwards, triangleheight(#_1, #_2, point) <= 0);
			);


			or(resultForwards, resultBackwards);
		);
				


		// *************************************************************************************************
		// Computes the angle at q from p to r. The result lies in [0, 2 * pi].
		// *************************************************************************************************
		computeangle(p, q, r) := (
		    regional(x, y, s, w);

		     x = p - q;
		     y = r - q;
		     s = (x * y) / (abs(x) * abs(y));
		     s = if(s < -1, -1, if(s > 1, 1, s));
		     w = arccos(s) + 0;

		     if(perp(x) * y >= 0, w, 2*pi - w);
		);


		// *************************************************************************************************
		// The sign of a number x.
		// *************************************************************************************************
		sign(x) := if(abs(x) == 0, 0, x / abs(x));


		// *************************************************************************************************
		// The faculty of the positive number n.
		// *************************************************************************************************
		faculty(n) := if(n <= 0, 1, n * faculty(n - 1));

		





		poissonDiskSampling(rect, d, numberOfPoints, searchThreshold) := (
			regional(oldPoints, hSize, vSize, result, searching, i, j, candidate, candidateValid, rangeA, rangeB);

			hSize = ceil(rect.w / d);
			vSize = ceil(rect.h / d);
			oldPoints = const((hSize + 1) * (vSize + 1), 0);

			result = [];
			
			searching = true;
			candidateValid = true;
			numberOfSearches = 0;

			while(and(searching, length(result) < numberOfPoints),
				candidate = [random() * rect.w, random() * rect.h];
				i = floor(candidate.x / d);
				j = floor(candidate.y / d);

				rangeA = max(i - 1, 0)..min(i + 1, hSize);
				rangeB = max(j - 1, 0)..min(j + 1, vSize);
				
				

				forall(rangeA, a, forall(rangeB, b,
					
					if(oldPoints_(a * vSize + b + 1) != 0, 
						candidateValid = and(candidateValid,
							dist(candidate, oldPoints_(a * vSize + b + 1)) > d;
						);	
					);
				));

				if(candidateValid,
					oldPoints_(i * vSize + j + 1) = candidate;
					result = result :> candidate;
					numberOfSearches = 0;
				, // else //		
					numberOfSearches = numberOfSearches + 1;
					if(numberOfSearches > searchThreshold,
						searching = false;	
					);
					candidateValid = true;
				);
			);

			apply(result, # + rect.xy);
		);
		poissonDiskSampling(rect, d, numberOfPoints) := poissonDiskSampling(rect, d, numberOfPoints, 32);







		// *************************************************************************************************
		// Checks, whether two line segments intersect.
		// *************************************************************************************************
		intersect(a, b) := (
		      area(a_1, a_2, b_1) * area(a_1, a_2, b_2) < 0
		    & area(b_1, b_2, a_1) * area(b_1, b_2, a_2) < 0
		);


		
		// *************************************************************************************************
		// Computes the signed distance of x to the line a-b.
		// *************************************************************************************************
		triangleheight(a, b, x) := if(a ~= b, dist(x, a), det(a :> 1, b :> 1, x :> 1) / dist(a, b));




		// *************************************************************************************************
		// Computes logarithms to p digits after the decimal point.
		// *************************************************************************************************
		arbiLog(x, b, p) := 10^(-p) * round(10^p * log(x) / log(b));
		log2(x, p)       := arbiLog(x, 2, p);
		log2(x)          := log2(x, 6);



		// *************************************************************************************************
		// The Kaylee way to scale probabilities.
		// *************************************************************************************************
		fuzzyScale(x, c) := c * x / ((c - 1) * x + 1);
		fuzzyS(x, c, lambda) := (x < lambda, lambda * fuzzyScale(x / lambda, 1 / c),  lambda + (1 - lambda) * fuzzyScale((x - lambda) / (1 - lambda), c));

		fuzzyVectorScale(x, c) := apply(x, fuzzyScale(#, c));


		
		// *************************************************************************************************
		// Computes the fractional part of a float.
		// *************************************************************************************************
		residual(x) := x - floor(x);
		fract(x)    := residual(x);



		// *************************************************************************************************
		// Computes a random point on the unit circle.
		// *************************************************************************************************
		randomDirection() := (
			regional(alpha);
			
			alpha = 2*pi*random();
			[sin(alpha), cos(alpha)];
		);


		// *************************************************************************************************
		// Returns the middle of three numbers.
		// Call it as mid(x, min_value, max_value) to clamp x in the interval [min_value, max_value].
		// *************************************************************************************************
		mid(a, b, c) := [a,b,c]--[min([a,b,c]),max([a,b,c])];



		// *************************************************************************************************
		// Checks, whether the distance of a point is less than eps to a segment (given by two points).
		// *************************************************************************************************
		closeToSegment(p, seg, eps) := (
		  (abs(triangleheight(seg_1, seg_2, p)) < eps)
		  & ((seg_2 - seg_1) * (p - seg_1) >= -0.5 * eps * abs(seg_2 - seg_1) * abs(p - seg_1))
		  & ((seg_1 - seg_2) * (p - seg_2) >= -0.5 * eps * abs(seg_1 - seg_2) * abs(p - seg_2));
		);




		// *************************************************************************************************
		// Linear interpolation between x and y.
		// *************************************************************************************************
		lerp(x, y, t) := t * y + (1 - t) * x;
		inverseLerp(x, y, p) := if(dist(y, x) != 0, (p - x) / (y - x), 0);
		// Lerp relative to t in interval [a, b].
		lerp(x, y, t, a, b) := lerp(x, y, inverseLerp(a, b, t));

		slerp(u, v, t) := (
			regional(angle);

			angle = arccos(u * v);

			if(angle ~= 0,
				u;
			, // else //
				(sin((1 - t) * angle) * u + sin(t * angle) * v) / sin(angle);
			);
		);
		inverseSlerp(u, v, w) := (
			regional(result);

			result = abs(min(inverseLerp(arctan2(u.x, u.y), arctan2(v.x, v.y), arctan2(w.x, w.y)), inverseLerp(arctan2(u.x, u.y) + 2*pi, arctan2(v.x, v.y), arctan2(w.x, w.y))));

			if(result > 1, 2 - result, result);

		);
		
		eerp(x, y, t) := x^(1 - t) * y^t;
		inverseEerp(x, y, p) := inverseLerp(log(x), log(y), log(p));





		// *************************************************************************************************
		// Computes signed area of simple polygons (i.e. not self intersecting.)
		// *************************************************************************************************
		areaOfPolygon(points) := (
		  0.5 * abs(sum(apply(1..(length(points) - 1),
		    (points_#).x * (points_(# + 1)).y - (points_(# + 1)).x * (points_#).y;
		  )));
		);




		lissajoue(a, b, d, t) := [sin(a * t + d), cos(b * t)];





        
		// *************************************************************************************************
		// Given a list of points, this returns the box circumscribing the points.
		// *************************************************************************************************
		box(list) := (
		    regional(bl, diag);

		    bl   = min(list);
		    diag = max(list) - bl;
		    expandrect(bl, diag.x, diag.y);
		);


        























		// *************************************************************************************************
		// Creates a rectangle (as a list of its vertices) of width w and height h with at position pos.
		// The value c gives the position in relation to the rectangle according to the number pad on a key-
		// board: E.g. 2 means the position is in the center of the bottom side, 5 means the center and 7
		// means the top left corner.
		// No c results in positioning the rectangle with pos in its lower left corner. I.e. c = 1 by
		// default.
		// *************************************************************************************************
		expandrect(pos, c, w, h) := (
		    regional(d, e, shift);

		    d     = 0.5 * [w, h];
		    e     = (d_1, -d_2);
		    shift = -compass()_c;
		    shift = (0.5 * w * shift.x, 0.5 * h * shift.y);
		    apply([-d, e, d, -e], pos + # + shift); //LU, RU, RO, LO
		);
		expandrect(pos, w, h) := expandrect(pos, 1, w, h);
		// Uses rect object; see below.
		expandrect(r) := expandrect(r.xy, r.c, r.w, r.h);

		compass() := -apply(directproduct([1, 0, -1], [1, 0, -1]), reverse(#));
		
		
		
		// ************************************************************************************************
		// Draws a rectangle with rounded corners.
		// ************************************************************************************************
		roundedrectangle(tl, w, h, r) := roundedrectangle(tl, tl + [w,-h], r);
		roundedrectangle(tl, br, r) := (
			regional(tr, bl);
			tr = [br.x, tl.y];
			bl = [tl.x, br.y];
			r = min([r, |tl.x-br.x|/2, |tl.y-br.y|/2]);
			//rounded corners
			circle(tl.xy + [r,-r], r)
				++ circle(bl.xy + [r,r], r)
				++ circle(br.xy + [-r,r], r)
				++ circle(tr.xy + [-r,-r], r)
			//rectangle
				++ polygon([tl.xy + [r,0], tr.xy + [-r,0], br.xy + [-r,0], bl.xy + [r,0]])
				++ polygon([tl.xy + [0,-r], tr.xy + [0,-r], br.xy + [0,r], bl.xy + [0,r]]);
		);


		
		
		// ************************************************************************************************
		// Draws and handles button. They have to be a JSON with the following keys and value-types:
		// button = {
		//   "position":   (2D vector),
		//   "size":       (2D vector),
		//   "label":      (String),
		//   "textSize":   (float),
		//   "colors":     (array with 3 colour vectors),
		//   "corner":     (float),
		//   "pressed":    (bool),
		//   "fontFamily": (String)
		// };
		// ************************************************************************************************
		drawButton(button) := (
		    if(button.pressed,
		        fill(roundedrectangle(button.position + 0.5 * (-button.size.x, button.size.y), button.size.x, button.size.y, button.corner), color -> (button.colors)_2);
		        fill(roundedrectangle(button.position + 0.5 * (-button.size.x, button.size.y) + (0, -0.2), button.size.x, button.size.y, button.corner), color -> (button.colors)_1);
		            drawtext(button.position + (0, -0.5 * button.textSize / 35) + (0, -0.2), button.label, align->"mid", size->button.textSize, color->(1, 1, 1), bold->true, family->button.fontFamily);
		    , // else //
		        fill(roundedrectangle(button.position + 0.5 * (-button.size.x, button.size.y) + (0, -0.2), button.size.x, button.size.y, button.corner), color -> (button.colors)_3);
		        fill(roundedrectangle(button.position + 0.5 * (-button.size.x, button.size.y), button.size_1, button.size_2, button.corner), color -> (button.colors)_2);
		        drawtext(button.position + (0, -0.5 * button.textSize / 35), button.label, align->"mid", size->button.textSize, color->(1, 1, 1), bold->true, family->button.fontFamily);
		    );
		);


		// Boilerplate code for button functionality. Call via
		// if(mouseInButton(button),
		// 	SOME CODE
		// );
		// in the mousedownscript and mouseupscript. The property button.pressed has to be set/updated manually; allowing for both switch- and toggle-buttons.
		mouseInButton(button) := (
		  dist(mouse().x, button.position.x) < 0.5 * button.size.x
		& dist(mouse().y, button.position.y) < 0.5 * button.size.y;
		);




		// ************************************************************************************************
		// Draws and handles toggles. They have to be a JSON with the following keys and value-types:
		// toggle = {
		//   "position":   (2D vector),
		//   "radius":     (float),
		//   "lineSize":   (float),
		//   "label":      (String),
		//   "textSize":   (float),
		//   "color":      (3D Vector),
		//   "pressed":    (bool),
		//   "fontFamily": (String)
		// };
		// ************************************************************************************************
		drawToggle(toggle) := (
			fillcircle(toggle.position, toggle.radius, size->toggle.lineSize, color->(1,1,1));
			if(toggle.pressed,
				drawcircle(toggle.position, toggle.radius, size->5 * toggle.lineSize, color->toggle.color);
			, // else //
				drawcircle(toggle.position, toggle.radius, size->toggle.lineSize, color->(0,0,0));
			);
			drawtext(toggle.position + [0, -0.015 * toggle.textSize], toggle.label, size->toggle.textSize, align->"mid", family->toggle.fontFamily);
		  );

		  catchToggle(toggle) := if(dist(mouse().xy, toggle.position) < toggle.radius, toggle.pressed = !toggle.pressed);
  



		// *************************************************************************************************
		// Draws text with border.
		// *************************************************************************************************
		drawWithBorder(pos, txt, size, align, color, bordercolor, bordersize, family) := (
			forall(bordersize * apply(1..8, [sin(2 * pi * #/ 8), cos(2 * pi * #/ 8)]), o,
				   drawtext(pos, txt, color -> bordercolor, offset -> o, size -> size, align -> align, family -> family);
				  );
			drawtext(pos, txt, color -> color, size -> size, align -> align, family -> family);
		  );
		  drawWithBorder(pos, txt, size, align, color, bordercolor, bordersize) := (
			forall(bordersize * apply(1..8, [sin(2 * pi * #/ 8), cos(2 * pi * #/ 8)]), o,
				   drawtext(pos, txt, color -> bordercolor, offset -> o, size -> size, align -> align);
				  );
			drawtext(pos, txt, color -> color, size -> size, align -> align);
		  );



		// *************************************************************************************************
		// Draws text over multiple lines with given line length.
		// *************************************************************************************************
		drawWrappedText(pos, txt, lineLength, lineHeight, size, align, color, family) := (
			regional(lines);

			lines = splitString(txt, lineLength);
			forall(1..length(lines),
				drawtext(pos + [0, -(# - 1) * lineHeight], lines_#, size->size, align->align, color->color, family->family);
			);
		);
		splitString(string, subLength) := (
			regional(totalSplit, result, candidate);

			totalSplit = tokenize(string, " ");

			result = [""];


			while(length(totalSplit) > 0,
				candidate = result_(-1) + totalSplit_1;
				
				if(length(candidate) <= subLength,
					result_(-1) = candidate + " ";
				, // else //
					result = result :> totalSplit_1 + " ";
				);

				totalSplit = bite(totalSplit);
			);

			result;
		);






		// *************************************************************************************************
		// Handles typing effect for LaTeX formulas.
		// *************************************************************************************************
		convertTexDelimiters(string, buffer) := (
			regional(startIndex, endIndex, innerCount, foundEnd);
		  
			startIndex = indexof(string, texDelimitersEBOW_1);
			if(startIndex == 0,
			  string;
			, // else //
			  endIndex = startIndex;
			  innerCount = 0;
			  foundEnd = false;
			  while(not(foundEnd),
				endIndex = endIndex + 1;
				if(endIndex == length(string),
				  foundEnd = true;
				, // else //
				  if(string_endIndex == texDelimitersEBOW_1, innerCount = innerCount + 1);
				  if(string_endIndex == texDelimitersEBOW_2, 
					innerCount = innerCount - 1;
					if(innerCount == -1,
					  foundEnd = true;            
					);
				  );
				)
			  );
		  
			  convertTexDelimiters(
				  sum(string_(1..startIndex - 1)) 
				+ if(buffer > 0, "", "\\phantom{")
				+ sum(string_(startIndex + 1 .. endIndex - 1))
				+ if(buffer > 0, "", "}") 
				+ sum(string_(endIndex + 1 .. length(string)))
			  , buffer - 1);
			  
		  
			);
		  );
		  parseTex(string) := apply(0..frequency(tokenize(string, ""), texDelimitersEBOW_1), convertTexDelimiters(string, #));

		  

		// ************************************************************************************************
		// Draws tex formula. Needs to have
		// toggle = {
		//   "pos":        (2D vector),
		//   "offset":     (2D vector)
		//   "pyramid":    (Array of strings), 
		//   "size":       (float),
		//   "size":       (float),
		//   "alpha":      (float),
		//   "color":      (3D Vector),
		//   "align":      (String)
		// };
		// ************************************************************************************************
		drawTexObject(obj) := drawtext(obj.pos + obj.offset, obj.outputString, size->obj.size, color->obj.color, align->obj.align, alpha->obj.alpha);

		createTexObject(pos, inputString, size) := (
			regional(pyramid);
			
			pyramid = parseTex(inputString);

			{
				"pos": pos,
				"offset": [0,0],
				"inputString": inputString,
				"pyramid": pyramid,
				"outputString": pyramid_1,
				"finalString": pyramid_(-1),
				"progress": 0,
				"size": size,
				"align": "left",
				"alpha": 1,
				"color": (0,0,0)
			};
		);
		typeTex(obj, track) := (
			obj.progress = clamp(track.progress, 0, 1);
			obj.outputString = obj.pyramid_(round(lerp(1, length(obj.pyramid), obj.progress)));
		);






		// *************************************************************************************************
		// Creates and handles rectangle objects.
		// *************************************************************************************************
		rect(x, y, c, w, h) := {"x": x, "y": y, "w": w, "h": h, "c": c, "xy": [x, y]};
		rect(x, y, w, h) := rect(x, y, 1, w, h);
		rect(pos, w, h)  := rect(pos.x, pos.y, w, h);
		drawRect(rect, size, color, alpha) := drawpoly(expandrect(rect), size -> size, color -> color, alpha -> alpha);
		fillRect(rect, color, alpha) := fillpoly(expandrect(rect), color -> color, alpha -> alpha);

		updateRectPos(rect, x, y) := (
			rect.x = x;
			rect.y = y;
			rect.xy = [x,y];
		);
		updateRectPos(rect, pos) := (
			rect.x = pos.x;
			rect.y = pos.y;
			rect.xy = pos;
		);

		rectContainsPoint(rect, point) := (
		  regional(expanded);

		  expanded = expandrect(rect);

		    point.x > (expanded_1).x
		  & point.x < (expanded_2).x
		  & point.y > (expanded_1).y
		  & point.y < (expanded_3).y
		);












		// *************************************************************************************************
		// Creates and handles slider UI element. Has to be a JSON object with the following keys and value-types.
		// slider = {
		//   "position":    (2D vector),
		//   "length":      (float),
		//   "size":        (float),
		//   "vertical":    (bool),
		//   "color":       (color vector),
		//   "startLabel":  (string),
		//   "endLabel":    (string),
		//   "labelSize":   (float),
		//   "value":       (float),
		//   "bulbSize":    (float),
		//   "fontFamily":  (string)
		// };
		// *************************************************************************************************
		sliderEnds(slider) := [slider.position, slider.position + if(slider.vertical, [0, -slider.length], [slider.length, 0])];

		drawSlider(slider) := (
		  regional(endPoints, startOffest, endOffset);

		  endPoints = sliderEnds(slider);

		  draw(endPoints, size -> slider.size, color -> slider.color);
		  fillcircle(lerp(endPoints_1, endPoints_2, slider.value), slider.bulbSize, color -> slider.color);
		  fillcircle(lerp(endPoints_1, endPoints_2, slider.value), 0.7 * slider.bulbSize, color -> (1,1,1));

		  startOffset = if(slider.vertical,
		    [0, 1.2 * slider.bulbSize + 0.2];
		  , // else //
		    [-1.2 * slider.bulbSize, -0.015 * slider.labelSize];
		  );
		  endOffset = if(slider.vertical,
		    [0, -1.2 * slider.bulbSize - 0.05 * slider.labelSize];
		  , // else //
		    [1.2 * slider.bulbSize, -0.015 * slider.labelSize];
		  );
		  drawtext(endPoints_1 + startOffset, slider.startLabel, align -> if(slider.vertical, "mid", "right"), size -> slider.labelSize, family->slider.fontfamily);
		  drawtext(endPoints_2 + endOffset,   slider.endLabel,   align -> if(slider.vertical, "mid", "left"),   size -> slider.labelSize, family->slider.fontfamily);
		);

		catchSlider(slider) := (
		  regional(endPoints);

		  endPoints = sliderEnds(slider);

		  if(closeToSegment(mouse().xy, endPoints, slider.bulbSize),
		    slider.value = if(slider.vertical,
		      inverseLerp((endPoints_1).y, (endPoints_2).y, mouse().y);
		    , // else //
		      inverseLerp((endPoints_1).x, (endPoints_2).x, mouse().x);
		    );
		  );
		);

		// *************************************************************************************************
		// Creates and handles selector UI element. Has to be a JSON object with the following keys and value-types.
		// selector = {
		//   "position":    (2D vector),
		//   "length":      (float),
		//   "size":        (float),
		//   "vertical":    (bool),
		//   "color":       (color vector),
		//   "textSize":    (float),
		//   "content":			(array),
		//   "index":       (int),
		//   "bulbSize":    (float),
		//   "fontFamily":  (String)
		// };
		// *************************************************************************************************
		drawSelector(selector) := (
			regional(endPoints);

		  endPoints = sliderEnds(selector);

		  draw(endPoints, size -> selector.size, color -> selector.color);
		  fillcircle(lerp(endPoints_1, endPoints_2, selector.index, 1, length(selector.content)), selector.bulbSize, color -> selector.color);
		  fillcircle(lerp(endPoints_1, endPoints_2, selector.index, 1, length(selector.content)), 0.7 * selector.bulbSize, color -> (1,1,1));

			forall(1..length(selector.content),
		    drawwithborder(lerp(endPoints_1, endPoints_2, #, 1, length(selector.content)) + (0, -0.015 * selector.textSize), selector.content_#, selector.textSize, "mid", [0,0,0], [1,1,1], 1.5, selector.fontFamily);
		  );

		);

		catchSelector(selector) := (
			regional(endPoints);

			endPoints = sliderEnds(selector);

			if(closeToSegment(mouse().xy, endPoints, selector.bulbSize),
		  	forall(1..length(selector.content),
			    if(dist(mouse().xy, lerp(endPoints_1, endPoints_2, #, 1, length(selector.content))) < selector.bulbSize,
						selector.index = #;
					);
			  );
		  );
		);



        ///////////////////////////////////////////////////////////////////////////////////////////////////
	// This library is intended to be used with CindyGL! But it doesnt have to...
	// It also needs egdodMath to work!
	///////////////////////////////////////////////////////////////////////////////////////////////////

		// *************************************************************************************************
		// Standard smmothstep function.
		// *************************************************************************************************
		smoothstep(x) := x * x * (3 - 2 * x);



		// *************************************************************************************************
		// Computes a random value in the interval [0,1] at a x,y position.
		// *************************************************************************************************
		randomValue(pos) := fract(sin(pos * [12.9898, 78.233]) * 43758.5453123);

		randomGradient(pos) := [2 * fract(sin(pos * (127.1,311.7)) * 43758.5453) - 1, 
			                    2 * fract(sin(pos * (269.5,183.3)) * 43758.5453) - 1 ];
		randomPoint(pos) := [fract(sin(pos * (127.1,311.7)) * 43758.5453), 
							 fract(sin(pos * (269.5,183.3)) * 43758.5453) ];

		// *************************************************************************************************
		// Gives random gradient noise based on a point in the plane.
		// *************************************************************************************************
		perlinNoise(coords) := (
			regional(iPoint, fPoint);
			
			iPoint = [floor(coords.x), floor(coords.y)];
			fPoint = [fract(coords.x), fract(coords.y)];

			0.5 * lerp(
					lerp(
						randomGradient(iPoint) * (fPoint), 
						randomGradient(iPoint + [1,0]) * (fPoint - [1,0]), 
					smoothstep(fPoint.x)),
					lerp(
						randomGradient(iPoint + [0,1]) * (fPoint - [0,1]),
						randomGradient(iPoint + [1,1]) * (fPoint - [1,1]),
					smoothstep(fPoint.x)),
				smoothstep(fPoint.y)) + 0.5;
		);
		perlinNoiseOctaves(coords) := (
			sum(apply(0..2, pow(0.5, #) * perlinNoise(pow(2, #) * coords))) / 1.75;
		);

		// *************************************************************************************************
		// Gives noise based on Voronoi decomposition.
		// *************************************************************************************************
		voronoiNoise(coords) := (
			regional(iPoint, mDist, neighbours, currDist);

			iPoint = [floor(coords.x), floor(coords.y)];
			fPoint = [fract(coords.x), fract(coords.y)];

			neighbours = directproduct(-1..1, -1..1);

			mDist = 1.5;

			forall(neighbours,
				currDist = dist(coords, iPoint + # + randomPoint(iPoint + #));	
				if(currDist < mDist,
					mDist = currDist;
				);
			);

			mDist;

		);
















        
		// ************************************************************************************************
		// Gets the ALICE and Toolbox colours.
		// ************************************************************************************************
		toolColorTable = dict(blue      -> [(2, 88,155),    (0, 70,124),      (0, 50, 89),     (0, 27, 48),     (0,  5,  9)]   / 255,
							  red       -> [(255,193,148),  (247,156, 90),    (227,114, 34),   (172, 76,  8),   (114, 47,  0)]    / 255,
							  green     -> [(237,247, 81),  (215,228, 15),    (162,173,  0),   (119,127,  0),   (61, 65,  0)]   / 255,
							  violet    -> [(219,90, 159),  (191,66, 144),    (168,37,  132),  (143,43,  132),  (116, 45,  133)]   / 255,
							  grey      -> [(224,224,224),  (172,172,172),    (128,128,128),   (78,78,78),      (35,35,35)]   / 255
		);
		toolColor(name, bright) := get(toolColorTable, name)_(3 - bright);
		toolColor(name)         := get(toolColorTable, name)_3;

		aliceColorTable = dict(blue      -> [(82,136,188),  (48,110,171),  (12,90,166),   (9,69,128),   (5,54,100)]   / 255,
							   orange    -> [(255,185,99),  (255,166,57),  (255,140,0),   (197,108,0),  (155,85,0)]   / 255,
							   violet    -> [(179,109,212), (158,62,204),  (145,19,204),  (123,5,178),  (95,3,138)]   / 255,
							   turquoise -> [(90,177,171),  (24,156,146),  (0,128,119),   (0,95,88),    (0,60,56)]    / 255,
							   teal      -> [(90,177,171),  (24,156,146),  (0,128,119),   (0,95,88),    (0,60,56)]    / 255,
							   red       -> [(253,105,109), (236,62,66),   (215,19,24),   (173,6,10),   (133,0,3)]    / 255,
							   green     -> [(131,222,92),  (99,204,53),   (67,186,16),   (48,150,5),   (34,116,0)]   / 255,
							   minion    -> [(255,242,153), (255,237,114), (245,224,80),  (221,198,44), (170,151,22)] / 255,
							   brown     -> [(186,104,43),  (152,74,16),   (117,50,0),    (80,34,0),    (42,18,0)]    / 255,
							   pink      -> [(255,246,255), (253,181,253), (248,131,248), (240,90,240), (225,52,225)] / 255,
							   grey      -> [(224,224,224), (172,172,172), (128,128,128), (78,78,78),   (36,35,35)]   / 255
		);
		aliceColor(name, bright) := get(aliceColorTable, name)_(3 - bright);
		aliceColor(name)         := get(aliceColorTable, name)_3;




		aliceColor = {
			"blue": (12,90,166) / 255,
			"orange": (255,140,0) / 255,
			"violet": (145,19,204) / 255,
			"turquoise": (0,128,119) / 255,
			"teal": (0,128,119) / 255,
			"red": (215,19,24) / 255,
			"green": (67,186,16) / 255,
			"minion": (245,224,80) / 255,
			"brown": (117,50,0) / 255,
			"pink": (248,131,248) / 255,
			"grey": (128,128,12) / 255
		};


		toolColor = {
			"blue": (0, 50, 89) / 255,
			"red": (227,114, 34) / 255, 
			"green": (162,173,  0) / 255, 
			"violet": (168,37,  132) / 255,
			"grey": (128,128,128) / 255 
		};

		sapColor = {
			"black":       hex2rgb("000000"),
			"darkGrey":    hex2rgb("7f7f7f"),
			"grey":        hex2rgb("aaaaaa"),
			"lightGrey":   hex2rgb("d4d4d4"),
			"white":       hex2rgb("ffffff"),
			"orange":      hex2rgb("ff8c00"),
			"lightOrange": hex2rgb("ffe54f"),
			"lightGreen":  hex2rgb("b0eb44"),
			"green":       hex2rgb("43ba10"),
			"lightBlue":   hex2rgb("39c0d5"),
			"blue":        hex2rgb("0c5aa6"),
			"background":  hex2rgb("1a2345"),
			"violet":      hex2rgb("9e3ecc"),
			"lightViolet": hex2rgb("ee9cff"),
			"lightRed":    hex2rgb("fd6886"),
			"red":         hex2rgb("d71318")
		};

    






		canvasRect = rect(canvasCorners.bl, canvasWidth, canvasHeight);












	
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
