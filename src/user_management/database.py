from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy.ext.declarative import declarative_base

engine = create_engine('sqlite:////Library/WebServer/Documents/BigSuR/databases/misc.db')
db_session = scoped_session(sessionmaker(autocommit=False,
                                         autoflush=False,
                                         bind=engine))
Base = declarative_base()
Base.query = db_session.query_property()

def init_db():
    # import all modules here that might define models so that
    # they will be registered properly on the metadata.  Otherwise
    # you will have to import them first before calling init_db()
    import user_management.models
    from user_management.models import User
    import uuid

    Base.metadata.create_all(bind=engine)

    # Ensure existing users have an fs_uniquifier (required by Flask-Security >=4.0.0)
    sess = None
    try:
        sess = db_session()
        need_update = False
        for u in sess.query(User).all():
            if not getattr(u, 'fs_uniquifier', None):
                setattr(u, 'fs_uniquifier', str(uuid.uuid4()))
                need_update = True
        if need_update:
            sess.commit()
    except Exception:
        # If the DB isn't initialized yet or there are issues, ignore here
        if sess is not None:
            sess.rollback()
    finally:
        if sess is not None:
            sess.close()