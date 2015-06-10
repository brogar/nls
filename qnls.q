\d .nls
loaded: 0b;

nls: `nls 2:(`q_nls;6);

defaultOpts: (`maxfev;`nprint) ! (2000; 0);

loaded:1b;
\d .


nlsfit:{[df;model;options]
	df: () xkey df;
	options: .nls.defaultOpts,options;
	options: key[options] ! `float$ value[options];
	model[`xvars]: (),model`xvars;
    df: @[df; model[`xvars],model`yvar; {`float$x}];
	
	Y: (0N;1) # ?[df;();();model`yvar];
	fn: raze string model`function;
	X: flip value ?[df;();(); model[`xvars]!model`xvars];
	Th: (0N;1) # model`initial;
	W: (count[Y];1) # `float$1.0;
	if[`weights in key model; W:(0N;1) # `float$ ?[df;();();model`weights]];
	
	call: .[ .nls.nls ; (Y;fn;X;Th;W;options); {'x}];
	call: `info`coef`resid`dof!call;
	:call;
    };

