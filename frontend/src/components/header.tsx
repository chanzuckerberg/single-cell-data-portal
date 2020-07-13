import { Link } from "gatsby";
import React, { FC } from "react";
import { Box, Flex, Image } from "theme-ui";
import Authenticate from "./authenticate";

const HomepageLink = () => {
  return (
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
  );
};

const Container: FC = props => {
  return (
    <Flex
      {...props}
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
    />
  );
};

const Header: FC = () => {
  return (
    <Box
      sx={{
        background: `#f5f5f5`,
        marginBottom: "0px",
      }}
    >
      <Container>
        <HomepageLink />
        <Authenticate />
      </Container>
    </Box>
  );
};

export default Header;
