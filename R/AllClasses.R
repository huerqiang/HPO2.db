HPODb <- setRefClass("HPODb", contains="AnnotationDb")
dbQuery <- getFromNamespace("dbQuery", "AnnotationDbi")


setMethod("keys", "HPODb",
    function(x, keytype, ...){
        if(missing(keytype)) keytype <- "hpoid"
        toTable(HPOTERM)[, keytype]
    }
)



setMethod("keytypes", "HPODb",
    function(x) {
        c("hpoid", "term")
    }

)

setMethod("select", "HPODb",
    function(x, keys, columns, keytype = "hpoid", ...){
        keytype <- match.arg(keytype, c("hpoid","term"))
        ## 获取hpoid
        strKeys <- paste0("\"", keys, "\"", collapse = ",")
        if (keytype == "term") {
            sql_key <- paste("SELECT hpoid FROM do_term WHERE term in (",
                strKeys, ")")
            hpoids <- dbQuery(dbconn(x), sql_key)[, 1]
            strKeys <- paste0("\"", hpoids, "\"", collapse = ",")
        }
        columns <- unique(c("hpoid", columns))

        sqls <- paste("SELECT ", paste(columns, collapse = ","),
            " FROM do_term")
        columns2 <- setdiff(columns, c("hpoid", "term"))
        for (col in columns2) {
            leftJoin <- paste0("LEFT JOIN  ", paste0("do_",col,
                " USING (hpoid)"))
            sqls <- c(sqls, leftJoin)
        }
        sqls <- c(sqls, paste0("WHERE do_term.hpoid in (", strKeys, ")"))
        sqls <- paste(sqls, collapse = " ")
        res <- dbQuery(dbconn(x), sqls)
        res
    }
)



setMethod("columns", "HPODb",
    function(x) {
        c("hpoid","term", "alias", "synonym", "parent", "children",
            "ancestor", "offspring")
    }
)
