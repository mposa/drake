Vreset = V2;
close(2)
figure(2)
hold on
%%
for i=11:100,
  Vreset = resetTest(x,s,Vreset,r,reset_constraint,V);
  
  if mod(i,2) == 0,
    contourSpotless(Vreset,x(1),x(3),[-.5 .5],[-1 1],[t;x([2;4;5])],zeros(model.num_states,1),1,{'r'});
  else
    contourSpotless(Vreset,x(1),x(3),[-.5 .5],[-1 1],[t;x([2;4;5])],zeros(model.num_states,1),1,{'b'});
  end
  
  
end