search_for_files <- function(file)
{
    if (file.exists(paste0("gene/", file, ".txt.gz")) &
        file.exists(paste0("variant/", file, ".txt.singleAssoc.txt.gz"))) {
        gene_file <- paste0("gene/", file, ".txt.gz")
        variant_file <- paste0("variant/", file, ".txt.singleAssoc.txt.gz")
        gz <- TRUE
    } else if (file.exists(paste0("gene/", file, ".txt")) &
               file.exists(paste0("variant/", file, ".txt.singleAssoc.txt"))) {
        gene_file <- paste0("gene/", file, ".txt")
        variant_file <- paste0("variant/", file, ".txt.singleAssoc.txt")
        gz <- FALSE
    } else if (file.exists(paste0("gene/", file, "_F.txt.gz")) &
               file.exists(paste0("variant/", file, "_F.txt.singleAssoc.txt.gz"))) {
        gene_file <- paste0("gene/", file, "_F.txt.gz")
        variant_file <- paste0("variant/", file, "_F.txt.singleAssoc.txt.gz")
        gz <- TRUE
    } else if (file.exists(paste0("gene/", file, "_F.txt")) &
               file.exists(paste0("variant/", file, "_F.txt.singleAssoc.txt"))) {
        gene_file <- paste0("gene/", file, "_F.txt")
        variant_file <- paste0("variant/", file, "_F.txt.singleAssoc.txt")
        gz <- FALSE
    } else {
        return(NULL)
    }
    return(list(gene_file=gene_file, variant_file=variant_file, gz=gz))
}

add_N <- function(file_info, dt_gene)
{
    # Read in the variant file to extract the sample size
    dt_variant <- fread(cmd = ifelse(file_info$gz,
        paste0("gzcat ", file_info$variant_file),
        paste0("cat ", file_info$variant_file)))

    if ("N" %in% colnames(dt_variant)) {
        N_eff <- mean(dt_variant$N)
        dt_gene$N_eff <- N_eff
        binary <- FALSE
    } else {
        N_case <- mean(dt_variant$N_case)
        N_control <- mean(dt_variant$N_ctrl)
        N_eff <- (4 * N_case * N_control) / (N_case + N_control)
        dt_gene$N_eff <- N_eff
        dt_gene$N_case <- N_case
        dt_gene$N_control <- N_control
        binary <- TRUE
    }
    return(list(dt_gene = dt_gene, binary=binary))
}

cauchy_combination <- function(p_values, weights=NULL)
{
	is.zero <- sum(p_values == 0) >= 1
	is.one <- sum(p_values > (1 - 1e-14)) >= 1

	if (is.zero) {
		return(0)
	}

	if (is.one) {
		p <- min(p_values) * length(p_values)
		if (p > 1) {
			return(-Inf)
		} else {
			return(qcauchy(p, lower.tail=FALSE))
		}
		return(min(1, (min(p_values)) * (length(p_values))))
	}

	if (is.null(weights)) {
		weights <- rep(1 / length(p_values), length(p_values))
	} else if (length(weights) != length(p_values)) {
		stop("The length of weights should be the same as that of the p-values!")
	} else if (sum(weights < 0) > 0) {
		stop("All the weights must be positive!")
	} else {
		weights <- weights / sum(weights)
	}

	is_small <- (p_values < 1e-16)
	if (sum(is_small) == 0) {
		cct_stat <- sum(weights * tan((0.5 - p_values) * pi))
	} else {
		cct_stat <- sum((weights[is_small] / p_values[is_small]) / pi)
		cct_stat <- cct_stat + 
			sum(weights[!is_small] * tan((0.5 - p_values[!is_small]) * pi))
	}

	return(cct_stat)
}

weighted_fisher <- function(p_values, weights=NULL, two_tail=FALSE, input_beta=NULL)
{
	if (is.null(weights)) {
		weights <- rep(1, length(p_values))
	}

	idx.na <- which(is.na(p_values))
	
	if (length(idx.na) > 0) {
		p_values <- p_values[-idx.na]
		weights <- weights[-idx.na]
		if (two_tail) {
			input_beta <- input_beta[-idx.na]
		}
	}

	NP <- length(p_values)
	NS <- length(weights)
	if (NP != NS) { stop("The length of p and weights vector must be identical.") }

	N <- NS
	Ntotal <- sum(weights)
	ratio <- weights / Ntotal
	Ns <- N * ratio
	G <- c()

	if (!two_tail) {
		for (i in 1:length(p_values)) {
		  G <- append(G, qgamma(p = p_values[i], shape = Ns[i], scale = 2, lower.tail=FALSE))
		}
		Gsum <- sum(G)
		resultP <- pgamma(q = Gsum, shape = N, scale = 2, lower.tail = FALSE)
	} else {
		p1 <- p2 <- p_values
		idx_pos <- which(input_beta > 0)
		idx_neg <- which(input_beta < 0)
		
		# Positive direction
		G <- c()
		p1[idx_pos] <- p_values[idx_pos] / 2
		p1[idx_neg] <- 1 - p_values[idx_neg] / 2

		for (i in 1:length(p1)) {
		  G <- append(G, qgamma(p = p1[i], shape = Ns[i], scale = 2, lower.tail = FALSE))
		}
		Gsum <- sum(G)
		resultP1 <- pgamma(q = Gsum, shape = N, scale = 2, lower.tail = FALSE)
		
		# Negative direction
		G <- c()
		p2[idx_pos] <- 1 - p_values[idx_pos] / 2
		p2[idx_neg] <- p_values[idx_neg] / 2

		for (i in 1:length(p2)) {
		  G <- append(G, qgamma(p = p2[i], shape = Ns[i], scale = 2, lower.tail = FALSE))
		}
		Gsum <- sum(G)
		resultP2 <- pgamma(q = Gsum, shape = N, scale = 2, lower.tail = FALSE)
		resultP <- 2 * min(resultP1, resultP2)
		if (resultP > 1.0) {
			resultP <- 1.0
		}
	}
	return(min(1, resultP))
}

run_weighted_fisher <- function(
	grouped_dt, n_eff_name, input_pvalues, output_meta_pvalue,
	two_tail = FALSE, input_beta = NULL)
{
	if (two_tail) {
		result <- grouped_dt %>% 
		summarise(
			"{output_meta_pvalue}" := weighted_fisher(
				.data[[input_pvalues]],
				weights = if (is.null(n_eff_name)) { NULL } else { .data[[n_eff_name]] },
				two_tail = two_tail, input_beta=.data[[input_beta]]
			)
		)
	} else {
		result <- grouped_dt %>% 
		summarise(
			"{output_meta_pvalue}" := weighted_fisher(
				.data[[input_pvalues]],
				weights = if (is.null(n_eff_name)) { NULL } else { .data[[n_eff_name]] }
			)
		)
	}
	return(result)
}

het_test <- function(input_beta, weights=NULL)
{

	beta_meta=sum_betas/sum_weights
	het_p=het_test(effs_size_org, weights, beta_meta)

    n_studies <- length(input_beta)
    effect_size_deviations  <- weights * (input_beta - effect_size_meta)^2
    return(pchisq(sum(effect_size_deviations, n_studies-1, lower.tail=FALSE)))
}

run_heterogeneity <- function(
	grouped_dt, n_eff_name, input_beta, output_meta_beta)
{
	grouped_dt <- grouped_dt %>%
	mutate(
		weights = sqrt(.data[[n_eff_name]]),
		beta = .data[[input_beta]]
	)
	summary_dt <- grouped_dt %>% 
	summarise(
		sum_weights = sum(weights),
		sum_betas = sum(weights * beta)
		) %>% mutate("{output_meta_beta}" := sum_betas/sum_weights)
	return(
		merge(grouped_dt, summary_dt) %>% 
		mutate(
			deviation = weights * (beta - .data[[output_meta_beta]])^2
		) %>% 
		summarise(sum_deviation = sum(deviation), n_studies=n()) %>% 
		mutate(p_het = pchisq(sum_deviation, n_studies-1))
	)
}

run_fisher <- function(
	grouped_dt, chi2_stat_name, input_pvalues, output_meta_pvalue)
{
	return(
		grouped_dt %>% 
		summarise(df = 2*n(),
			"{chi2_stat_name}" := -2*sum(log(.data[[input_pvalues]]))) %>%
		mutate("{output_meta_pvalue}" := pchisq(.data[[chi2_stat_name]],
			df=df, lower.tail=FALSE)) %>% select(-df)
	)
}

run_cauchy <- function(
	grouped_dt, n_eff_name, Cauchy_stat_name, input_pvalues, output_meta_pvalue) 
{
	return(
		grouped_dt %>% 
		summarise(
			"{Cauchy_stat_name}" := cauchy_combination(
				.data[[input_pvalues]],
				weights=sqrt(.data[[n_eff_name]])
			),
			number_of_pvals := n()
		) %>% mutate(
		  	"{output_meta_pvalue}" := ifelse(
				.data[[Cauchy_stat_name]] > 1e+15,
				(1 / .data[[Cauchy_stat_name]]) / pi,
				pcauchy(.data[[Cauchy_stat_name]], lower.tail=FALSE)
			)
		) %>% mutate(
			"{output_meta_pvalue}" := ifelse(
				.data[[output_meta_pvalue]] > (1 - 1e-10),
				(1 - 1/number_of_pvals), .data[[output_meta_pvalue]]
			)
		)
	)
}

run_inv_var <- function(
	grouped_dt, beta_meta_name, se_meta_name,
	input_beta, input_SE, output_meta_pvalue
) {
	return(
		grouped_dt %>% 
		mutate(weight = 1/(.data[[input_SE]]**2)) %>%
		mutate(effs_inv_var = .data[[input_beta]] * weight) %>%
		summarise(
			"{beta_meta_name}" := sum(effs_inv_var) / sum(weight),
			"{se_meta_name}" := sqrt(1/sum(weight)),
			"{output_meta_pvalue}" := 2 * dnorm(
				abs(.data[[beta_meta_name]] / .data[[se_meta_name]]))
		)
	)
}

run_stouffer <- function(
	grouped_dt, n_eff_name, weighted_Z_name,
	input_pvalues, output_meta_pvalue,
	two_tail = FALSE, input_beta = NULL
) {
	if (two_tail) {
		grouped_dt <- grouped_dt %>% 
			mutate("{input_pvalues}" := .data[[input_pvalues]]/2)
	} else {
		input_beta <- "beta_dummy"
		grouped_dt <- grouped_dt %>% mutate("{input_beta}" := 1)
	}
	
	result <- grouped_dt %>%
		mutate(
			weighted_Z_numerator = (
				sqrt(.data[[n_eff_name]]) * 
				(-qnorm(.data[[input_pvalues]])) * 
				sign(.data[[input_beta]])
			)
		)

	if (two_tail) {
		result <- result %>%
		summarise(
			"{weighted_Z_name}" := sum(weighted_Z_numerator) / 
				sqrt(sum(.data[[n_eff_name]])),
			"{output_meta_pvalue}" := 2 * pnorm(abs(.data[[weighted_Z_name]]), lower.tail=FALSE)
		)
	} else {
		result <- result %>%
		summarise(
			"{weighted_Z_name}" := sum(weighted_Z_numerator) / 
				sqrt(sum(.data[[n_eff_name]])),
			"{output_meta_pvalue}" := pnorm(.data[[weighted_Z_name]], lower.tail=FALSE)
		)
	}
	return(result)
}