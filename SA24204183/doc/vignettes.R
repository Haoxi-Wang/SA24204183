## -----------------------------------------------------------------------------
library(SA24204183)

# Example datasets
train_data <- rnorm(200)  # Example training data (normal distribution)
test_data <- rnorm(50)    # Example test data (normal distribution)

# Call KL_for_sales function with default parameters
result <- KL_for_sales(i = 1, train = train_data, test = test_data)

# Display results
result

## -----------------------------------------------------------------------------
# Load the package (replace 'YourPackageName' with your actual package name)
library(SA24204183)

# Set the number of iterations and thinning steps
N <- 1000  # Total number of samples
thin <- 10  # Thinning step

max(N, thin)

# Generate samples using the gibbsC function
samples <- gibbsC(N = N, thin = thin)

# Display the first few samples of x and y
head(samples)

