########## Nutrigonometry IV: Thales' theorem #################
## Scientific Reports
## Juliano Morimoto
#* Please cite the original paper if using the approach or the functions

# Adapted by AM Senior 2023

#*30/03/2023*#

## Functions ####

# Functions
calculate_hypothenuse <- function(x, y){
  sqrt(x^2 + y^2)
}

calculate_largeTriangle <- function(data, intake = TRUE, intake_coord = NULL){
  if(isTRUE(intake) & !is.null(intake_coord)){
    intake_coorddt <- intake_coord
    h <- calculate_hypothenuse(intake_coorddt[[1]], intake_coorddt[[2]])
    s <- calculate_hypothenuse(data[[1]], data[[2]])
    t <- calculate_hypothenuse((data[[1]] - intake_coorddt[[1]]), 
                               (data[[2]] - intake_coorddt[[2]]))
    
    out_tri <- cbind(h, s, t)
    colnames(out_tri) <- c("diameter_h", "s", "t")
    return(as.data.frame(out_tri))
  } else {
    stop('Intake coordinates must be given')
  }
}



## Calculate Beta

calculate_Beta <- function(dt){
  ## estimating angles
  output <- c()
  for(i in 1:nrow(dt)){
    # using Cosine rule #
    beta <- acos((dt$s[i]^2 + dt$t[i]^2 - (dt$diameter_h[i]^2))/(2*dt$s[i]*dt$t[i]))*57.2958
    output[i] <- beta 
  }
  return(output)
}


calculate_Thales <- function(data, intake_coord, ...){
  ## calculate variables 
  tri_out <- calculate_largeTriangle(data = data, intake = TRUE, intake_coord = intake_coord)
  ## Calculate angle
  out_beta <- calculate_Beta(tri_out)
  return(out_beta)
}

# Function to double check calcualtion of beta angle
my_beta<-function(X, Y, IT_X, IT_Y){
	
	# Written out in full	
	# acos((X^2 + Y^2 + (X - IT_X)^2 + (Y - IT_Y)^2 - (IT_X^2 + IT_Y^2)) / (2 * sqrt(X^2 + Y^2) * sqrt((X - IT_X)^2 + (Y - IT_Y)^2))) * 57.2958
	
	# Simplifies to
	out<-acos((Y * (Y - IT_Y) - IT_X * X + X^2) / (sqrt(X^2 + Y^2) * sqrt((X - IT_X)^2 + (Y - IT_Y)^2))) * 57.2958
	return(out)
}

# Function to get SE for beta by delta method
my_SE<-function(int_x, int_y, diet, IT_x, IT_y){
	
	# Load the msm package
	require(msm)
	
	# Formula for beta as above
	form<-sprintf("~ acos((x2 * (x2 - %.0f) - %.0f * x1 + x1^2) / (sqrt(x1^2 + x2^2) * sqrt((x1 - %.0f)^2 + (x2 - %.0f)^2))) * 57.2958", IT_y, IT_x, IT_x, IT_y)
		
	# Aggregate in to a data set
	aggregated<-data.frame(x=int_x, y=int_y, diet=diet)
	diets<-unique(aggregated$diet)
	SE<-NULL
	
	# Apply delta method to each diet
	for(i in 1:length(diets)){
		subset<-aggregated[which(aggregated$diet==diets[i]),]
		mod<-lm(cbind(x, y) ~ 1, data=subset)
		mu<-mod$coef
		V<-vcov(mod)
		se<-deltamethod(g=as.formula(form), mean=mu, cov=V)
		SE<-c(SE, se)
	}
	
	# Parcel up and return
	return(data.frame(diet=diets, se=SE))
	
}

# Function to get and mean se with bootstrapping
boot_mu<-function(dat, indicies){
	d<-dat[indicies]
	return(c(mean(d), sd(d)/sqrt(length(d)-1)))
}


# Caitlin's functions for food intake


##### My CDM Model

CDM<-function(foods, IT){
	
	# Make sure the food_comps are a proportion
	f<-function(x){
		x/sum(x)
	}
	food_comp<-t(apply(foods, 1, f))

	# Calculate the amounts of each nutrient eaten on a food to reach the IT
	f<-function(x){
		V<-c(-IT %*% x) / sum(x^2)
		abs(V) * x
	}
	intakes<-(t(apply(food_comp, 1, f)))
	return(intakes)	

}




