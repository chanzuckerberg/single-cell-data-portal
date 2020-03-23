import { Link } from "gatsby"
import React from "react"
import { Flex, Box, Image } from "theme-ui"

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
          maxWidth: "1100px",
          padding: [3],
          paddingLeft: 5,
        }}
      >
        <Box sx={{ width: 112 }}>
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
      </Flex>
    </Box>
  )
}

export default Header
