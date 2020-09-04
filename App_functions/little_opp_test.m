function out = little_opp_test (in)
test_nr = numel(in);
for ii = 1:test_nr-1
   if in(ii) ~= in(ii+1)
       continue
   else
       in(ii) = NaN;
   end
       
end
in_log = ~isnan(in);

out = in(in_log);

end
