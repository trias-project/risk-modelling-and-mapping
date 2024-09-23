
conns <- showConnections(all = TRUE)

# Loop through and close each open connection
for (i in seq_len(nrow(conns))) {
  con_id <- as.integer(rownames(conns)[i])
  if (con_id > 2) {  # Exclude standard connections 0 (stdin), 1 (stdout), 2 (stderr)
    close(getConnection(con_id))
  }
}


# Optionally, check if all connections are closed
showConnections(all = TRUE)
