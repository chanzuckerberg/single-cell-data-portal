import React, { FC } from "react";
import { Link } from "gatsby";
import { Box, Flex, Image } from "theme-ui";
import Authenticate from "./authenticate";

const Header: FC = () => {
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
          maxWidth: 1,
          paddingTop: 0,
          paddingBottom: 0,
          paddingLeft: 5,
          paddingRight: 5,
        }}
      >
        <Link
          to="/"
          style={{
            color: `black`,
            textDecoration: `none`,
          }}
        >
          <Image
            src="https://chanzuckerberg.com/wp-content/themes/czi/img/logo-minified.svg"
            width="40px"
            height="40px"
            marginTop="5px"
          />
        </Link>
        <Authenticate />
      </Flex>
    </Box>
  );
};

export default Header;
