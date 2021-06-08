message('Loading oncoKB API functions')
library(httr)
########## load oncoKB API functions #############

#' @name get_oncokb_gene_entry_url
#' @title get_oncokb_gene_entry_url
#' @description
#'
#' Get the url for the OncoKB entry for an actionable alteration
#'
#' @param oncokb.response the JSON response, which is an output of the httr::GET query
#' @return The URL if exists, otherwise NA
#' @export
#' @author Alon Shaiber
get_oncokb_gene_entry_url = function(oncokb.response){
    if (!inherits(oncokb.response, 'list')){
        oncokb.response = list(oncokb.response)
    }
    urls = sapply(oncokb.response, function(response){
        baseurl = 'https://www.oncokb.org/gene'
        if (!inherits(response, 'response')){
            stop('You must provide an response object or a list of response objects.')
        }
        oncokb.content = httr::content(response)
        if (oncokb.content$geneExist == FALSE){
            return(NA)
        }
        gene = strsplit(oncokb.content$geneSummary, ',')[[1]][1]
        if (length(gene) == 0 | !is.character(gene)){
            return(NA)
        }
        if (!(length(oncokb.content$treatments) == 0)){
        # Notice: I am using the treatments object to get the AA mutation
        # in the future we can switch to reading it from our VCF file
            alterations = oncokb.content$treatments[[1]]$alterations
            if (is.null(alterations)){
                return(NA)
            }
            if (!inherits(alterations, 'list')){
                return(NA)
            }
            alteration = alterations[[1]]
            if (is.character(alteration) & length(alteration) > 0){
                library(RCurl)
                # if a URL exists for the specific alteration then let's go there
                url_alteration = paste0(baseurl, '/', gene, '/', alteration, '/')
                if (url.exists(url_alteration)){
                    return(url_alteration)
                }
            }
        } else {
            # If there was no URL for the specific alteration then we will try to go to the gene
            url_gene = paste0(baseurl, '/', gene, '/')
            if (url.exists(url_gene)){
                return(url_gene)
            }
        }
        return(NA)
    })
    return(urls)
}


#' @name get_oncokb_response
#' @title get_oncokb_response
#' @description
#'
#' Using httr to query genomic changes on OncoKB in accordance with the web API specification in https://www.oncokb.org/swagger-ui/index.html .
#'
#' @param variants.dt data.table with the variants to query. As a minimum, must contain the following columns: seqnames, start, end, REF, ALT.
#' @param oncokb.token a token to use when querying OncoKB. If you don't have a token, you must first obtain one from: https://www.oncokb.org/apiAccess
#' @return list of response object. Response objects can be parsed using httr (for example: httr::content).
#' @export
#' @author Alon Shaiber
get_oncokb_response = function(variants.dt, oncokb.token,
                                  oncokb.url = 'https://www.oncokb.org/api/v1/annotate',
                                  reference = 'GRCh37'){
    message('Querying oncoKB')
    refdict = list('hg19' = 'GRCh37', 'hg38' = 'GRCh38')
    if (reference %in% names(refdict)){
        reference = refdict[[reference]]
    }
    if (!(reference %in% refdict)){
        stop('An invalid reference was provided. The only valid options are: ', paste(unique(c(names(refdict), refdict)), collapse = ', '))
    }
    if (!('start' %in% names(variants.dt))){
        if (!('pos' %in% names(variants.dt))){
            stop('variants.dt must either include a "pos" column or "seqnames", "start, and "end" columns')
        }
        tryCatch(
        {
            pos.dt = gr2dt(parse.gr(variants.dt$pos))
            variants.dt$seqnames = pos.dt$seqnames
            variants.dt$start = pos.dt$start
            variants.dt$end = pos.dt$end
        },
        error = function(e) print("The input variants.dt has a 'pos' column that is not valid. 'pos' must be parsable by parse.gr"))
    }
    all.res = lapply(
        1:nrow(variants.dt),
        function(i){
            res = GET(
                url = paste0(
                    oncokb.url, "/mutations/byGenomicChange?",
                    "genomicLocation=", variants.dt[i, as.character(seqnames)], "%2C",
                    variants.dt[i, start], "%2C",
                    variants.dt[i, end], "%2C",
                    variants.dt[i, REF], "%2C", variants.dt[i, ALT],
                    "&referenceGenome=", reference),
                add_headers(authorization = paste("Bearer", oncokb.token)),
                add_headers(accept = "application/json")
            )
        })
    return(all.res)
}

#' @name get_oncokb_annotations
#' @title get_oncokb_annotations
#' @description
#'
#' Get a data.table containing specific fields from a list of OncoKB responses
#'
#' @param oncokb.response either a single response object or a list of response objects
#' @param fields the fields to extract from the OncoKB response
#' @return data.table with the values of the fields in each response
#' @export
#' @author Alon Shaiber
get_oncokb_annotations = function(oncokb.response,
                                  fields = c('highestSensitiveLevel',
                                             'highestResistanceLevel',
                                             'highestDiagnosticImplicationLevel',
                                             'highestPrognosticImplicationLevel')){
    message('Pulling the oncoKB annotations from oncoKB responses')
    if (!inherits(oncokb.response, 'list')){
        oncokb.response = list(oncokb.response)
    }
    annotations = lapply(oncokb.response, function(response){
        if (!inherits(response, 'response')){
            stop('You must provide an response object or a list of response objects.')
        }
        error_code_received = FALSE
        oncokb_errors = NULL
        status = status_code(response)
        status_type = (status%/%100) * 100
        if (status_type != 200){
            error_code_received = TRUE
            reason = http_status(status)$reason
            oncokb_errors = sprintf("%s (HTTP %d).", reason, status)
            vdt = NULL
        } else {
            oncokb.content = httr::content(response)
            invalid.fields = setdiff(fields, names(oncokb.content))
            if (length(invalid.fields) > 0){
                stop('The following fields are not valid and are missing from the oncoKB response: ', paste(invalid.fields, collaps = ', '))
            }
            vdt = setNames(data.table(matrix(ncol = length(fields), nrow = 1)), fields)
            for (field in fields){
                value = oncokb.content[[field]]
                if (length(value) > 1){
                    stop('The following oncoKB field is not valid: ', field, ' since it contains more than one item.')
                }
                if (!is.null(value)){
                    vdt[, (field) := value]
                }
            }
        }
        return(list(vdt, error_code_received, oncokb_errors))
        })
    print(annotations)
    print('hi')
    print(lapply(annotations, '[', 2))
    print('HI')
    print(lapply(annotations, '[', 3))
    error_code_received = any(unlist(sapply(annotations, '[', 2)))
    if (error_code_received){
        oncokb_errors = paste(unique(sapply(annotations, '[', 3)), collapse = ', ')
        message('The following errors were received from oncoKB: ', oncokb_errors)
        message('Skipping oncoKB annotations.')
        return(NULL)
    }
    annotations = rbindlist(lapply(annotations, '[', 1), fill = TRUE)
    return(annotations)
}

#' @name test_oncokb_api
#' @title test_oncokb_api
#' @description
#'
#' A function to test our internal API with OncoKB
#'
#' @author Alon Shaiber
test_oncokb_api = function(oncokb.token){
    # this is an example actionable genomic alteration
    example1 = data.table(seqnames = '7', start = '140453136', end = '140453136', REF = 'A', ALT = 'T')
    # this is a gene that exists in the database, but a variant that does not
    example2 = data.table(seqnames = '2', start = '204732714', end = '204732714', REF = 'A', ALT = 'G')
    # this is a gene that is not in the database at all
    example3 = data.table(seqnames = '1', start = '69511', end = '69511', REF = 'A', ALT = 'G')
    som = rbind(example1, example2, example3)
    oncokb = get_oncokb_response(som, oncokb.token = oncokb.token)
    oncokb_annotations = get_oncokb_annotations(oncokb)
    message('Get the entry url for oncoKB')
    onco_kb_entry_url = get_oncokb_gene_entry_url(oncokb)
    if (is.na(oncokb_annotations[1, highestSensitiveLevel])) stop('Something went wrong. The query for a known actionable genomic alteration did not return the expected value.')
    message('OncoKB hg19 query was finished with success.')

    # example using GRCh38
    example_GRCh38 = data.table(seqnames = 'chr7', start = '140753336', end = '140753337', REF = 'A', ALT = 'T')
    oncokb38 = get_oncokb_response(example_GRCh38, reference = 'GRCh38', oncokb.token = oncokb.token)
    oncokb_annotations = get_oncokb_annotations(oncokb38)
    onco_kb_entry_url = get_oncokb_gene_entry_url(oncokb38)
    if (is.na(oncokb_annotations[1, highestSensitiveLevel])) stop('Something went wrong. The query for a known actionable genomic alteration did not return the expected value.')
    message('OncoKB hg38 query was finished with success.')
}
########## done loading oncoKB API functions #############
