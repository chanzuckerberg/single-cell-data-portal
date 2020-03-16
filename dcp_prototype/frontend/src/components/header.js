import { Link } from "gatsby"
import React from "react"
import Logo from "./logo1"
import { Box, Flex } from "theme-ui"
import Authenticate from "./authenticate"

const Header = ({ siteTitle }) => {
  return (
    <Box
      sx={{
        background: `#f5f5f5`,
        marginBottom: [4],
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
          <Link
            to="/"
            style={{
              color: `black`,
              textDecoration: `none`,
            }}
          >
            <Logo/>
          </Link>
        </Box>
        <Box>Data</Box>
        <Box>Guides</Box>
        <Box>Metadata</Box>
        <Box>Pipelines</Box>
        <Box>Analysis Tools</Box>
        <Box>Contribute</Box>
        <Box>APIs</Box>
        <Box>|</Box>
        <Authenticate/>
      </Flex>
    </Box>
  )
}

export default Header
