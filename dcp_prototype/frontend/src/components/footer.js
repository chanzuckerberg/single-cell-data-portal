import { Link } from "gatsby"
import React from "react"
import Logo from "./logo2"
import { Flex, Box } from "theme-ui"

const Header = ({ siteTitle }) => {
  return (
    <Box
      sx={{
        background: `black`,
        mt: [4],
      }}
    >
      <Flex
        sx={{
          alignItems: "center",
          justifyContent: "space-between",
          fontWeight: 700,
          margin: `0 auto`,
          maxWidth: 960,
          padding: [3],
        }}
      >
        <Box sx={{ width: 200 }}>
          <Link to="/">
            <Logo />
          </Link>
        </Box>
        <Box sx={{ color: "white" }}>Data</Box>
        <Box sx={{ color: "white" }}>Guides</Box>
        <Box sx={{ color: "white" }}>Metadata</Box>
        <Box sx={{ color: "white" }}>Pipelines</Box>
        <Box sx={{ color: "white" }}>Analysis Tools</Box>
        <Box sx={{ color: "white" }}>Contribute</Box>
        <Box sx={{ color: "white" }}>APIs</Box>
      </Flex>
    </Box>
  )
}

export default Header
