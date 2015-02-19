function [w0_tild2,lambda_tild2,phi_tild2]=est_correction(w0_tild,lambda_tild,phi_tild,alpha_tild2,beta_tild2,p)
    w0_tild2 = w0_tild-(alpha_tild2*beta_tild2)/p;
    lambda_tild2 = lambda_tild-0.25*alpha_tild2/p+0.25*log(1+(beta_tild2/p)^2);
    phi_tild2 = phi_tild+0.25*(alpha_tild2)^2*beta_tild2/p-0.5*atan(beta_tild2/p);
end