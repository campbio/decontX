library(decontX)

## Test internal error checking
mat <- matrix(seq(5), ncol = 10, nrow = 10)
label1 <- rep(seq(2), each = 5)
label2 <- as.factor(seq(100))
label3 <- as.factor(label1)
label4 <- label3
label4[seq(2)] <- 2
label5 <- as.factor(rep(seq(5), each = 2))

test_that(desc = "Testing .colSumByGroup", {
    expect_error(.Call("_colSumByGroup", mat, label1))
    expect_error(.Call("_colSumByGroup", mat, label2))
    res <- .Call("_colSumByGroup", mat, label3)
    expect_true(all(res == t(rowsum(t(mat), label3))))
    res <- .colSumByGroup(mat, label3, 2)
    expect_true(all(res == t(rowsum(t(mat), label3))))
})

storage.mode(mat) <- "numeric"

test_that(desc = "Testing .colSumByGroupNumeric", {
    expect_error(.Call("_colSumByGroup_numeric", mat, label1))
    expect_error(.Call("_colSumByGroup_numeric", mat, label2))
    res <- .Call("_colSumByGroup_numeric", mat, label3)
    expect_true(all(res == t(rowsum(t(mat), label3))))
    res <- .colSumByGroupNumeric(mat, label3, 2)
    expect_true(all(res == t(rowsum(t(mat), label3))))
})
