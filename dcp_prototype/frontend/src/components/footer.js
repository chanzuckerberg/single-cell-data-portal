import { Link } from "gatsby"
import React from "react"
import Logo from "./logo2"
import { Flex, Box } from "theme-ui"

const Header = ({ siteTitle }) => {
  return (
    <Box
      sx={{
        background: `black`,
        mt: "0px",
      }}
    >
      <Flex
        sx={{
          alignItems: "center",
          justifyContent: "space-between",
          fontWeight: 700,
          margin: `0px`,
          maxWidth: "100%",
          padding: [3],
          paddingLeft: 5,
        }}
      >
        <Box sx={{ width: 112 }}>
          <a href="https://www.humancellatlas.org/">
            <Logo />
          </a>
        </Box>
        <Box sx={{ color: "white",
          display: `none`,
        }}>Data</Box>
        <Box sx={{ color: "white",
          display: `none`,
        }}>Guides</Box>
        <Box sx={{ color: "white",
          display: `none`,
        }}>Metadata</Box>
        <Box sx={{ color: "white",
          display: `none`,
        }}>Pipelines</Box>
        <Box sx={{ color: "white",
          display: `none`,
        }}>Analysis Tools</Box>
        <Box sx={{ color: "white",
          display: `none`,
        }}>Contribute</Box>
        <Box sx={{ color: "white",
          display: `none`,
        }}>APIs</Box>
      </Flex>
    </Box>
  )
}

export default Header
