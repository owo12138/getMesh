#' 获取 PubMed MeSH 主题词和条目词
#'
#' 这个函数从 PubMed 获取指定关键词的 MeSH 主题词和条目词，并返回它们。
#'
#' @param csv_file 字符串，包含输入数据的 CSV 文件路径，CSV 文件应该包含关键词列。
#'
#' @return 返回一个列表，包含 MeSH 主题词和条目词的向量。
#'
#' @examples
#' # 假设有一个名为 "keywords.csv" 的文件，其中包含关键词
#' pubmed_mesh_online("path_to_your_csv_file.csv")
#'
#' @export
pubmed_mesh_online <- function(csv_file = system.file("mesh.csv",package = "getMesh")) {

  # 函数：抓取 MeSH 主题词和 Entry Terms
  get_mesh_details <- function(keyword) {
    base_url <- "https://www.ncbi.nlm.nih.gov"
    search_url <- paste0(base_url, "/mesh/?term=", URLencode(keyword))

    # 使用 rvest 抓取网页内容
    cat("Debug: Accessing URL:", search_url, "\n")  # 调试信息
    page <- tryCatch({
      read_html(search_url)
    }, error = function(e) {
      stop("Error: Unable to access the MeSH search page.")
    })

    # 将抓取到的整个页面内容转换为字符串并保存到一个列表中以进行调试
    debug_page_list <<- append(debug_page_list, list(as.character(page)))

    # 提取主题词（假设主题词在一个特定的 HTML 标签中）
    topic <- page %>%
      html_nodes(xpath = '//*[@id="maincontent"]/div/div[5]/div/h1') %>%
      html_text(trim = TRUE)

    # 定位到包含 "Entry Terms:" 的 <p> 标签
    entry_terms_section <- page %>%
      html_nodes(xpath = '//*[text()="Entry Terms:"]/following-sibling::ul[1]')

    # 提取 Entry Terms 部分的内容
    entry_terms <- character(0)  # 初始化空的字符向量
    if (length(entry_terms_section) > 0) {
      # 提取 <ul> 标签中的所有 <li> 元素并将每个条目作为单独元素存储
      entry_terms <- entry_terms_section %>%
        html_nodes(xpath = './/li') %>%
        html_text(trim = TRUE)  # 提取文本内容并去除多余空格
    }
    return(list(topic = topic, entry_terms = entry_terms))
  }

  # 初始化调试页面列表
  debug_page_list <<- list()

  # 读取 CSV 文件
  data <- tryCatch({
    read_csv(csv_file)
  }, error = function(e) {
    stop("Error: Unable to read CSV file. Please check the file path and format.")
  })

  # 提取 MeSH 和条目词函数（使用在线抓取）
  extract_mesh_and_entry <- function(topics) {
    mesh_list_local <- map(topics, get_mesh_details)
    mesh_list_local <- mesh_list_local[!sapply(mesh_list_local, is.null)]

    if (length(mesh_list_local) == 0) {
      return(list(topics, topics))  # 如果找不到 MeSH 词，返回原始关键词作为 MeSH 和条目词
    }

    mh <- unlist(map(mesh_list_local, "topic"))
    entries <- unlist(map(mesh_list_local, "entry_terms"))

    # 创建并返回 list
    lists <- list(mh, entries)

    return(lists)
  }


  # 创建搜索查询函数
  create_search_query <- function(df) {
    topics <- unlist(df[2:length(df)])
    topics <- unique(topics[!is.na(topics) & topics != ""])

    if (length(topics) == 0) {
      warning("Warning: No valid topics found in this row.")
      return(NULL)
    }

    lists <- extract_mesh_and_entry(topics)
    mh_query <- unique(lists[[1]][lists[[1]] != ""])
    entry_query <- unique(lists[[2]][lists[[2]] != ""])

    # 将原始关键词加入查询
    original_topics_query <- paste0("\"", topics, "\"")

    # 去除重复项并组合所有查询
    queries <- unique(c(mh_query, entry_query, original_topics_query))

    # 移除空值
    queries <- queries[queries != ""]

    if (length(queries) == 0) {
      return(NULL)  # 如果生成的查询为空，返回 NULL
    }

    final_query <- create_nested_or_query(queries)
    return(final_query)
  }

  # 生成查询语句
  queries <- apply(data, 1, function(row) {
    create_search_query(row)
  })

  # 移除空值
  queries <- queries[!sapply(queries, is.null)]

  if (length(queries) == 0) {
    stop("Error: No valid queries generated from the input data.")
  }

  # 合并所有查询为最终的 PubMed 查询
  final_query <- Reduce(function(x, y) {
    paste0("(", x, " AND ", y, ")")
  }, queries)

  # 打印最终查询
  cat(final_query, sep = "\n")
}

#' Title
#'
#' @param terms
#'
#' @return
#'
#' @examples
create_nested_or_query <- function(terms) {
  terms <- unique(terms[terms != ""])  # 去除空字符串和重复项

  # 如果只剩下一个term，直接返回
  if (length(terms) == 1) {
    return(paste0("(", terms, ")"))
  }

  # 如果有多个term，递归构建深度嵌套的OR结构
  create_nested <- function(terms) {
    if (length(terms) == 1) {
      return(terms)
    }
    mid <- floor(length(terms) / 2)
    left <- create_nested(terms[1:mid])
    right <- create_nested(terms[(mid+1):length(terms)])
    return(paste0("(", left, " OR ", right, ")"))
  }

  return(create_nested(terms))
}
