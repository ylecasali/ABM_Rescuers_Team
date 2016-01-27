[payoff]=eval_payoff_rescued(injured_found,injured_rescued,payoff,val_found,val_rescued,val_else)

if (injured_found)
payoff = payoff + val_found;
elseif (injured_rescued)
payoff = payoff + val_rescued;
else
payoff = payoff - val_else;
end


