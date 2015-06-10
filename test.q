@[system;"l qnls.q";{'x}];

.Q.fs[{`DNase insert flip `Run`conc`density!("IFF";",")0:x}]`:data/DNase

DNase1: select from DNase where Run=1;

modelFun:{[x;theta] 
	/ all arguments are presumed lists
	conc: x[0];
	Asym: theta[0];
	xmid: theta[1];
	scal: theta[2];
	denr: 1.0 + exp[(xmid - log[conc])%scal];
	ans: Asym%denr;
	:ans;
	};

model: (`yvar;`xvars;`function;`initial) ! (`density; `conc; `modelFun; 3.0 -1.0 10.0);

modelOpts:(`xtol; `nprint; `gtol; `ftol; `maxfev)!(1.0e-9; 1f; 1.0e-10; 1.0e-11; 20);

fit: nlsfit[DNase1; model; modelOpts];
