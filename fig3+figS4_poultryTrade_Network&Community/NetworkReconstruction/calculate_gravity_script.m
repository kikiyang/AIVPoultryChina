%Live poultry trade
G_lptUpdate = calculate_gravity(poultryPopulation97_14,humanPopulation97_14,distance,1,1,1);
%Poultry egg trade
CGe2=calculate_gravity(eggProduction97_14*10000,humanPopulation97_14,distance,1,1,1);
%Live poultry trade - Cambodia Gravity Model Parameter
G_hk200=calculate_gravity(poultryPopulation97_14,humanPopulation97_14,distance,0.5427,0.9339,167.4561);