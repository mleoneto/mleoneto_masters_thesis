function gx = gModels2L(x_t,P,u,in)

model = in;

switch model
    case 'bias'
        
        alpha = P(1);
        temp = exp(P(2));
        diff_ut = alpha;
  
    case 'basic'
        true_dice = u(1);
        other_dice = u(2);
        
        alpha = exp(P(1));
        delta = P(2);
        temp = exp(P(3));
        
        diff_ut = (alpha.*(other_dice - true_dice) + (delta.*other_dice));
        
    case 'fixed_cost'
        true_dice = u(1);
        other_dice = u(2);
        alpha = exp(P(1));
        delta = P(2);
        temp = exp(P(3));

        diff_ut = (alpha.*(other_dice - true_dice) + delta);
        
    case 'action_alpha'
        alpha = exp(P(1));
        delta = P(2);
        temp = exp(P(3));
        
        true_dice = u(3);
        other_dice = u(4);
        bad_group = u(2);
        if bad_group == 1        
            alpha_up = x_t(2);
        else
            alpha_up = x_t(1);            
        end
        
        diff_ut = ((alpha + alpha_up).*(other_dice - true_dice) + (delta.*other_dice));
        
    case 'action_delta'
        alpha = exp(P(1));
        delta = P(2);
        temp = exp(P(3));
        
        true_dice = u(3);
        other_dice = u(4);
        bad_group = u(2);
        if bad_group == 1        
            delta_up = x_t(2);
        else
            delta_up = x_t(1);            
        end        
        
        diff_ut = ((alpha.*(other_dice - true_dice)) + (delta + delta_up).*other_dice);
        
    case 'action_alpha_delta'
        alpha = exp(P(1));
        delta = P(2);
        temp = exp(P(3));
        
        true_dice = u(3);
        other_dice = u(4);
        bad_group = u(2); 
        
        if bad_group == 0 
            alpha_up = x_t(1);
            delta_up = x_t(3);
        else
            alpha_up = x_t(2);
            delta_up = x_t(4);            
        end           
        
        diff_ut = ((alpha + alpha_up).*(other_dice - true_dice) + ((delta + delta_up).*other_dice));
        
    case 'outcome_delta'
        alpha = P(1);
        delta = P(2);
        temp = exp(P(3));
        
        true_dice = u(5);
        other_dice  = u(6);
        bad_group = u(2);
        
        if bad_group == 1        
            delta_up = x_t(2);
        else
            delta_up = x_t(1);            
        end            
        
        diff_ut = ((alpha.*(other_dice - true_dice)) + (delta + delta_up).*other_dice);
                
    case 'outcome_alpha'
        alpha = P(1);
        delta = P(2);
        temp = exp(P(3));
        
        true_dice = u(5);
        other_dice  = u(6);
        bad_group = u(2);        
        
        if bad_group == 1        
            alpha_up = x_t(2);
        else
            alpha_up = x_t(1);            
        end
        
        diff_ut = ((alpha + alpha_up).*(other_dice - true_dice) + (delta.*other_dice));
                        
    case 'outcome_alpha_delta'
        alpha = P(1);
        delta = P(2);
        temp = exp(P(3));
        
        true_dice = u(5);
        other_dice  = u(6);
        bad_group = u(2);        
        
        if bad_group == 0
            alpha_up = x_t(1);
            delta_up = x_t(3);
        else
            alpha_up = x_t(2);
            delta_up = x_t(4);            
        end          
        
        diff_ut = ((alpha + alpha_up).*(other_dice - true_dice) + ((delta + delta_up).*other_dice));
    
    case 'conformity'
        alpha = exp(P(1));
        delta = P(2);
        temp = exp(P(3));
        gamma = P(4);
        true_dice = u(2);
        other_dice  = u(3);
        basel = u(1);
        
        if basel == 1
        diff_ut = (alpha.*(other_dice - true_dice) + (delta.*other_dice));
        else 
        diff_ut = (alpha.*(other_dice - true_dice) + (delta.*other_dice)) + gamma;
        end
        
            case 'conformity_sim'
        alpha = P(1);
        delta = P(2);
        temp = exp(P(3));
        gamma = P(4);
        true_dice = u(1);
        other_dice  = u(2);
        
        diff_ut = (alpha.*(other_dice - true_dice) + (delta.*other_dice)) + gamma;
    
    case 'conformity_action'
       
        alpha = exp(P(1));
        delta = P(2);
        temp = exp(P(3));
        gamma = P(4);    
        true_dice = u(3);
        other_dice  = u(4);
        
         basel = u(5);
        
        if basel == 1
        diff_ut = (alpha.*(other_dice - true_dice) + (delta.*other_dice));
        else 
        diff_ut = (alpha.*(other_dice - true_dice) + (delta.*other_dice)) + gamma*x_t;
        end
        
        
    case 'conformity_outcome'
        alpha = exp(P(1));
        delta = P(2);
        temp = exp(P(3));
        gamma = P(4);    
        true_dice = u(5);
        other_dice  = u(6);
        basel = u(7);
        
        if basel == 1
        diff_ut = (alpha.*(other_dice - true_dice) + (delta.*other_dice));
        else 
        diff_ut = (alpha.*(other_dice - true_dice) + (delta.*other_dice)) + gamma*x_t;
        end        
end

gx = VBA_sigmoid(temp*diff_ut);

end
