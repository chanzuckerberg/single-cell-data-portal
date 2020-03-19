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
        marginBottom: "0px",
      }}
    >
      <Flex
        sx={{
          height: "60px",
          alignItems: "center",
          justifyContent: "space-between",
          fontWeight: 700,
          margin: `0 auto`,
          maxWidth: "100%",
          paddingTop: 0,
          paddingBottom: 0,
          paddingLeft: 5,
          paddingRight: 5,


        }}
      >
        <Box sx={{ width: 175 }}>
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
        <Box
					style={{
						display: `none`,
					}}>Data</Box>
        <Box
					style={{
						display: `none`,
					}}>Guides</Box>
        <Box
					style={{
						display: `none`,
					}}>Metadata</Box>
        <Box
					style={{
						display: `none`,
					}}>Pipelines</Box>
        <Box
					style={{
						display: `none`,
					}}>Analysis Tools</Box>
        <Box
					style={{
						display: `none`,
					}}>Contribute</Box>
        <Box
					style={{
						display: `none`,
					}}>APIs</Box>
      </Flex>
    </Box>
  )
}

export default Header
