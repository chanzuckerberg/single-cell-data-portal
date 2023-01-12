import { Auth0ContextInterface, useAuth0 } from "@auth0/auth0-react";
import { AnchorButton, Button, MenuDivider } from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { Menu, MenuItem } from "czifui";
import { MouseEventHandler, useState } from "react";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import { AuthButtonWrapper } from "src/components/Header/style";
import CuratorAPIKeyGenerator from "../CuratorAPIKeyGenerator";
import { LogOutEmail, LogOutText } from "./style";

const AuthButtons = (): JSX.Element | null => {
  const hasAuth = get(FEATURES.CURATOR) === BOOLEAN.TRUE;

  const {
    error,
    isAuthenticated,
    isLoading,
    loginWithRedirect,
    logout,
    user: userInfo,
  } = useAuth0();

  if (userInfo && error) {
    // (thuang): Force refresh page to log user out
    window.location.reload();
  }

  if (!hasAuth || isLoading) return null;

  return (
    <AuthButtonWrapper>
      {isAuthenticated ? (
        <LoggedInButtons
          name={userInfo?.name}
          email={userInfo?.email}
          logout={logout}
        />
      ) : (
        <LoggedOutButtons handleLogin={loginWithRedirect} />
      )}
    </AuthButtonWrapper>
  );
};

function LoggedInButtons({
  name,
  email,
  logout,
}: {
  name?: string;
  email?: string;
  logout: Auth0ContextInterface["logout"];
}) {
  const authName = isEmail(name) ? "Account" : name;
  const [anchorEl, setAnchorEl] = useState<Element | null>(null);

  const handleClick: MouseEventHandler<HTMLElement> = (event) => {
    setAnchorEl(event.currentTarget);
  };

  const handleClose = () => {
    setAnchorEl(null);
  };

  return (
    <>
      <Content />
      <Button
        onClick={handleClick}
        minimal
        rightIcon={IconNames.CHEVRON_DOWN}
        text={authName}
      />
    </>
  );

  function Content() {
    const curatorAPIFeature = get(FEATURES.CURATOR);

    return (
      <Menu
        anchorEl={anchorEl}
        keepMounted
        open={Boolean(anchorEl)}
        onClose={handleClose}
      >
        {curatorAPIFeature && (
          <div>
            <CuratorAPIKeyGenerator />
            <MenuDivider />
          </div>
        )}
        <MenuItem data-testid="log-out" onClick={handleLogout}>
          <LogOutText>Logout</LogOutText>
          <LogOutEmail data-testid="user-email">{email}</LogOutEmail>
        </MenuItem>
      </Menu>
    );
  }

  function handleLogout() {
    logout();
  }
}

function LoggedOutButtons({
  handleLogin,
}: {
  handleLogin: Auth0ContextInterface["loginWithRedirect"];
}) {
  return (
    <AnchorButton
      onClick={handleLogin}
      minimal
      rightIcon={IconNames.CHEVRON_DOWN}
      text="Log In"
    />
  );
}

function isEmail(name?: string): boolean {
  if (!name) return false;

  return name.includes("@");
}

export default AuthButtons;
