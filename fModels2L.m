function fx = fModels2L(x_t,P,u,in)

switch in
    case 'action_alpha'
        
        c_alpha1 = x_t(1);
        c_alpha2 = x_t(2);        
        %Get values of u
        diff = u(1);
        bad_group = u(2);
        
        %Get values of P
        lr = P(1);
        
        if bad_group == 1
            fx(2) = c_alpha2 + lr.*diff;
            fx(1) = c_alpha1;
        elseif bad_group == 0
            fx(1) = c_alpha1 + lr.*diff;
            fx(2) = c_alpha2;
        end
        
    case 'action_delta'
        
        c_delta1 = x_t(1);
        c_delta2 = x_t(2);       
        %Get values of u
        diff = u(1);
        bad_group = u(2);
        
        %Get values of P
        lr = P(1);
        
        if bad_group == 1
            fx(2) = c_delta2 + lr.*diff;
            fx(1) = c_delta1;
        elseif bad_group == 0
            fx(1) = c_delta1 + lr.*diff;
            fx(2) = c_delta2;
        end
        
    case 'action_alpha_delta'
        
        %Get current delta
        c_alpha1 = x_t(1); %good group
        c_alpha2 = x_t(2); %bad group
        c_delta1 = x_t(3);
        c_delta2 = x_t(4);        
        %Get values of u
        diff = u(1);
        bad_group = u(2);
        
        %Get values of P
        lr_alpha = P(1);
        lr_delta = P(2);
        
        if bad_group == 1
            fx(2) = c_alpha2 + lr_alpha.*diff;
            fx(4) = c_delta2 + lr_delta.*diff;
            fx(1) = c_alpha1;
            fx(3) = c_delta1;
        elseif bad_group == 0
            fx(1) = c_alpha1 + lr_alpha.*diff;
            fx(3) = c_delta1 + lr_delta.*diff;
            fx(2) = c_alpha2;
            fx(4) = c_delta2;            
        end
        
    case 'outcome_alpha'
        %Get current delta
        c_alpha1 = x_t(1);
        c_alpha2 = x_t(2);   
        
        %Get values of u
        truthful =  u(1) ;
        liars = 10 - truthful;
        true_dice = u(2) ;
        other_dice = u(3) ;
        bad_group = u(4) ;
        
        %Get values of P
        lr = P(1);
        
        if bad_group == 1
            fx(2) = c_alpha2 + lr.*((other_dice.*liars)/10);
            fx(1) = c_alpha1;
        elseif bad_group == 0
            fx(1) = c_alpha1 + lr.*((other_dice.*liars)/10);
            fx(2) = c_alpha2;
        end
        
    case 'outcome_delta'
        %Get current delta
        c_delta1 = x_t(1);
        c_delta2 = x_t(2);
        
        %Get values of u
        truthful =  u(1) ;
        true_dice = u(2) ;
        other_dice = u(3) ;
        bad_group = u(4) ;
        
        %Get values of P
        lr = P(1);
        
        if bad_group == 1
            fx(2) = c_delta2 + lr.*((true_dice.*truthful)/10);
            fx(1) = c_delta1;
        elseif bad_group == 0
            fx(1) = c_delta1 + lr.*((true_dice.*truthful)/10);
            fx(2) = c_delta2;
        end  
        
    case 'outcome_alpha_delta'
        %Get current delta
        c_alpha1 = x_t(1);
        c_alpha2 = x_t(2);        
        c_delta1 = x_t(3);
        c_delta2 = x_t(4);        
       
        %Get values of u
        truthful =  u(1) ;
        liars = 10 - truthful;
        true_dice = u(2) ;
        other_dice = u(3) ;
        bad_group = u(4) ;
        
        %Get values of P
        lr_alpha = P(1);
        lr_delta = P(2);
        
        if bad_group == 1
            fx(2) = c_alpha2 + lr_alpha.*((other_dice.*liars)/10);
            fx(4) = c_delta2 + lr_delta.*((true_dice.*truthful)/10);
            fx(1) = c_alpha1;
            fx(3) = c_delta1;            
        elseif bad_group == 0
            fx(1) = c_alpha1 + lr_alpha.*((other_dice.*liars)/10);
            fx(3) = c_delta1 + lr_delta.*((true_dice.*truthful)/10);
            fx(2) = c_alpha2;
            fx(4) = c_delta2;             
        end
    
  
     case 'conformity_action'
        diff =  u(1) ;
       init = u(2) ; 
     lr = P(1);
        
           
        if init == 1
            fx = diff;
        elseif init == 0
            fx = x_t + lr*diff;
        end 
        
    case 'conformity_outcome'
        diff =  u(1) ;
        init = u(2) ;  
        true_dice = u(3);
        other_dice = u(4) ;
       
        lr = P(1);
        
        if init == 1
            fx = diff*other_dice;
        elseif init == 0
            fx = x_t + lr*diff*(other_dice - true_dice);     
        end
             
       
                        
end

end
