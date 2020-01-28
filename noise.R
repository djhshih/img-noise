library(ggplot2)
library(magick)
library(filenamer)

sigma_x <- 1;
rho <- sqrt(0.3);

# By definition,
# \rho^2 = \frac{E[ (X - \mu_X) (Y - \mu_Y) ]^2}{\sigma^2_X}{\sigma^2_Y}

# Model
# X ~ Normal(0, \sigma^2_X)
# Y \mid X ~ Normal(X, \epsilon)

# Therefore, the marginal distribution of Y is
# Y ~ Normal(0, \sigma^2_Y)
# where \sigma^2_Y = \sigma^2_X + \epsilon

# In the noiseless setting, \epsilon = 0 => \sigma^2_Y = \sigma^2_X, and
# \rho^2_0 = \frac{E[ (X - \mu_X) (Y - \mu_Y) ]^2}{\sigma^2_X}{\sigma^2_X} = 1
# Therefore,
# \sigma^2_X = \frac{E[ (X - \mu_X) (Y - \mu_Y) ]^2}{\sigma^2_X}

# Substituting back to original equation gives
# \rho^2 = \frac{\sigma^2_X}{\sigma^2_X + \epsilon^2}
# Rearranging gives
# \epsilon = \sqrt{ \frac{ \sigma^_X (1 - \rho^2) }{ \rho^2 } }

compute_epsilon <- function(sigma_x, rho) {
	sqrt( sigma_x^2 * (1 - rho^2) / rho^2 );
}

# By definition,
# CV_Y = \frac{ \sigma_{Y - X} }{ \mu_Y }
# where \sigma^_{Y - X} = Var[Y - X] = \epsilon^2
# Therefore,
# \epsilon = CV_Y \mu_Y = CV_Y \mu_X

compute_epsilon_cv <- function(mu_x, cv) {
	cv * mu_x
}

epsilon <- compute_epsilon(sigma_x, rho);



n <- 10000;

x <- rnorm(n, 0, sigma_x);
y <- rnorm(n, x, epsilon);

sd(x)
sd(y)
sqrt(sigma_x^2 + epsilon^2)

print(cor(x, y))

####

# source: https://pixabay.com
#in.fname <- "img/cat-on-tree-small.jpg";
#in.fname <- "img/dog-in-forest-small.jpg";
in.fname <- "img/animals-in-forest-small.jpg";

out.fname <- as.filename(in.fname) %>% set_fpath(path=c("img", "noisy"));

img <- image_read(in.fname)
plot(img)

r <- as.integer(img[[1]]);
sd(r)
hist(r)
epsilon <- compute_epsilon(sd(r), rho);

epsilon2 <- compute_epsilon_cv(mean(r), 0.8)

r.n <- array(rnorm(length(r), r, epsilon), dim=dim(r));

sd(r.n)
hist(r.n)

cor(r, r.n)^2
sd(r.n - r) / mean(r)
smoothScatter(r, r.n)

# pixel wrapping

r.n.w <- round(r.n) %% 256;
sd(r.n.w)
hist(r.n.w)

cor(r, r.n.w)
smoothScatter(r, r.n.w)

# pixel bounding

r.n.b <- round(r.n);
r.n.b[r.n < 0] <- 0;
r.n.b[r.n > 255] <- 255;

sd(r.n.b)
hist(r.n.b)

cor(r, r.n.b)
smoothScatter(r, r.n.b)

# pixel rescaling

r.n.s <- (r.n - min(r.n)) / (diff(range(r.n)));

sd(r.n.s)
hist(r.n.s)

cor(r, r.n)^2
cor(r, r.n.s)^2
smoothScatter(r, r.n.s)

img.n.w <- image_read(as.raster(r.n.w, max=255));
plot(img.n.w)

img.n.b <- image_read(as.raster(r.n.b, max=255));
plot(img.n.b)

img.n.s <- image_read(as.raster(r.n.s));
plot(img.n.s)

plot(image_convert(img, colorspace="gray"))
plot(image_convert(img.n.s, colorspace="gray"))


image_add_noise_and_rescale <- function(img, rho, cv) {
	r <- as.integer(img[[1]]);
	if (!missing(rho)) {
		epsilon <- compute_epsilon(sd(r), rho);
	} else if (!missing(cv)) {
		epsilon <- compute_epsilon_cv(mean(r), cv);
	} else {
		stop("One of rho or cv must be specified")
	}
	r.n <- array(rnorm(length(r), r, epsilon), dim=dim(r));
	r.n.s <- (r.n - min(r.n)) / (diff(range(r.n)));

	image_read(as.raster(r.n.s))
}

# Whole-genome sequencing in high-confidence regions yielded 100%
# reproducibility across technical replicates:
# https://www.pnas.org/content/113/42/11901

rho2s <- c(0.9, 0.6, 0.3, 0.1);
rhos <- sqrt(rho2s);

imgs.n <- lapply(rhos, function(rho) image_add_noise_and_rescale(img, rho=rho));

plot(img)
plot(imgs.n[[1]])
plot(imgs.n[[2]])
plot(imgs.n[[3]])
plot(imgs.n[[4]])

plot(image_convert(imgs.n[[3]], colorspace="gray"))
plot(image_convert(imgs.n[[4]], colorspace="gray"))

lapply(imgs.n, function(x) cor(as.integer(x[[1]])[, , 1:3], as.integer(img[[1]]))^2)

dens <- 300;

make_path(out.fname);
mapply(
	function(i, rho) {
		image_write(i, tag(out.fname, sprintf("rho2-%.1f", rho)), format=out.fname$ext, density=dens)
		j <- image_convert(i, colorspace="gray");
		image_write(j, tag(out.fname, c(sprintf("rho2-%.1f", rho), "gray")), format=out.fname$ext, density=dens)
	},
	imgs.n, rho2s
)

cvs <- c(0.2, 0.5, 0.8, 1.4);

# CVs of reference genes in TCGA RNA expression data
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2809-2/tables/1

imgs.n.cv <- lapply(cvs, function(cv) image_add_noise_and_rescale(img, cv=cv));

plot(img)
plot(imgs.n.cv[[1]])
plot(imgs.n.cv[[2]])
plot(imgs.n.cv[[3]])
plot(imgs.n.cv[[4]])

plot(image_convert(imgs.n.cv[[4]], colorspace="gray"))

mapply(
	function(i, cv) {
		image_write(i, tag(out.fname, sprintf("cv-%.1f", cv)), format=out.fname$ext, density=dens)
		j <- image_convert(i, colorspace="gray");
		image_write(j, tag(out.fname, c(sprintf("cv-%.1f", cv), "gray")), format=out.fname$ext, density=dens)
	},
	imgs.n.cv, cvs
)

