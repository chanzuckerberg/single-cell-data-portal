import {
  AnchorButton,
  Button,
  Intent,
  Menu,
  MenuItem,
  Popover,
} from "@blueprintjs/core";
import { IconNames } from "@blueprintjs/icons";
import { API } from "src/common/API";
import { get } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { BOOLEAN } from "src/common/localStorage/set";
import { useUserInfo } from "src/common/queries/auth";
import { API_URL } from "src/configs/configs";
import { ButtonWrapper, Initial } from "./style";

const AuthButtons = () => {
  const hasAuth = get(FEATURES.CURATOR) === BOOLEAN.TRUE;

  const { data: userInfo, isLoading, error } = useUserInfo(hasAuth);

  if (userInfo && error) {
    // (thuang): Force refresh page to log user out
    window.location.reload();
  }

  if (!hasAuth || isLoading) return null;

  return (
    <ButtonWrapper>
      {userInfo?.name ? (
        <LoggedInButtons name={userInfo.name} email={userInfo.email} />
      ) : (
        <LoggedOutButtons />
      )}
    </ButtonWrapper>
  );
};

function LoggedInButtons({ name, email }: { name?: string; email?: string }) {
  const Name = isEmail(name) ? <Initial>{name ? name[0] : ""}</Initial> : name;

  return (
    <>
      <Popover content={<Content />}>
        <Button
          large
          outlined
          minimal
          intent={Intent.PRIMARY}
          rightIcon={IconNames.CARET_DOWN}
        >
          {Name}
        </Button>
      </Popover>
    </>
  );

  function Content() {
    return (
      <Menu>
        <MenuItem data-testid="user-email" text={`Logged in as: ${email}`} />
        <MenuItem
          data-testid="log-out"
          text="Log Out"
          href={`${API_URL}${API.LOG_OUT}`}
          icon={IconNames.LOG_OUT}
        />
      </Menu>
    );
  }
}

function LoggedOutButtons() {
  return (
    <>
      <AnchorButton
        href={`${API_URL}${API.LOG_IN}`}
        outlined
        intent={Intent.PRIMARY}
      >
        Log In
      </AnchorButton>
    </>
  );
}

function isEmail(name?: string): boolean {
  if (!name) return false;

  return name.includes("@");
}

export default AuthButtons;
