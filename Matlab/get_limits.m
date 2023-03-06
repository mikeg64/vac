% Calculate fmin and fmax from function values
% Uses nfunc, funcs, autorange, and isfirstpict variables.

for ifunc=1:nfunc
   if autorange(ifunc)
      fun=trim(funcs(ifunc,:));
      get_func;
      f_max=max2(f1);
      f_min=min2(f1);
      if isfirstpict
         fmax(ifunc)=f_max;
         fmin(ifunc)=f_min;
      else
         if f_max>fmax(ifunc);fmax(ifunc)=f_max;end;
         if f_min<fmin(ifunc);fmin(ifunc)=f_min;end;
      end
   end
end

