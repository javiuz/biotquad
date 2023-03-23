function sigma=build_initial_sigma

global lambda mu

sigma=zeros(48,1);

val=lambda+mu;

sigma(2)=val;
sigma(3)=val;

sigma(6)=val;
sigma(7)=val;
sigma(10)=val;

sigma(11)=val;
sigma(14)=val;

sigma(15)=val;
sigma(18)=val;
sigma(19)=val;

sigma(21)=val;
sigma(24)=val;
sigma(25)=val;
sigma(28)=val;

sigma(29)=val;
sigma(31)=val;
sigma(34)=val;

sigma(35)=val;
sigma(38)=val;

sigma(39)=val;
sigma(42)=val;
sigma(44)=val;

sigma(45)=val;
sigma(48)=val;
end