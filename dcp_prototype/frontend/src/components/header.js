import { Link } from "gatsby"
import React from "react"
import { Box, Flex, Image } from "theme-ui"
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
          maxWidth: "1100px",
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
          <Image src="https://chanzuckerberg.com/wp-content/themes/czi/img/logo-minified.svg"
                 width="55px"
                 height="55px"/>
          </Link>
        </Box>
        <Authenticate/>
      </Flex>
    </Box>
  )
}

export default Header
