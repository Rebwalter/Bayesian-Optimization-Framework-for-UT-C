function Params = defineParameters_opt()

Params.Height_canyon.lb = 15;
Params.Height_canyon.ub = 24.9;

Params.Width_canyon.lb  = 15.0;
Params.Width_canyon.ub  = 30.0;

Params.Width_roof.lb = 8.0;
Params.Width_roof.ub = 25.0;

Params.fveg_G.lb = 0.0;
Params.fveg_G.ub = 1;

% Intrinsic radius scaling
Params.alpha.lb = 0.02;     % Rt/Wc
Params.alpha.ub = 0.24;     % temporary upper, will shrink dynamically

Params.albedo_w.lb = 0.25;
Params.albedo_w.ub = 0.9;

Params.lan_dry_W.lb = 0.1;
Params.lan_dry_W.ub = 3.8;

Params.cv_s_W.lb = 0.1*10^6; 
Params.cv_s_W.ub = 3.2*10^6;

end